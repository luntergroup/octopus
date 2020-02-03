// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "median_base_quality.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MedianBaseQuality::name_ = "BQ";

std::unique_ptr<Measure> MedianBaseQuality::do_clone() const
{
    return std::make_unique<MedianBaseQuality>(*this);
}

Measure::ResultType MedianBaseQuality::get_default_result() const
{
    return std::vector<boost::optional<int>> {};
}

namespace {

bool is_canonical(const VcfRecord::NucleotideSequence& allele) noexcept
{
    const static VcfRecord::NucleotideSequence deleted_allele {vcfspec::deletedBase};
    return !(allele == vcfspec::missingValue || allele == deleted_allele);
}

bool has_called_alt_allele(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    if (!call.has_genotypes()) return true;
    const auto& genotype = get_genotype(call, sample);
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                       [&] (const auto& allele) { return allele != call.ref() && is_canonical(allele); });
}

bool is_evaluable(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    return has_called_alt_allele(call, sample);
}

auto median_base_quality(const ReadRefSupportSet& reads, const Allele& allele)
{
    boost::optional<AlignedRead::BaseQuality> result {};
    if (!is_indel(allele)) {
        std::vector<AlignedRead::BaseQuality> base_qualities {};
        base_qualities.reserve(reads.size() * sequence_size(allele));
        for (const auto& read : reads) {
            if (overlaps(read.get(), allele)) {
                utils::append(copy_base_qualities(read, mapped_region(allele)), base_qualities);
            }
        }
        if (!base_qualities.empty()) {
            result = maths::median(base_qualities);
        }
    }
    return result;
}

} // namespace

Measure::ResultType MedianBaseQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    std::vector<boost::optional<int>> result {};
    result.reserve(call.num_alt());
    for (const auto& sample : samples) {
        boost::optional<int> sample_result {};
        if (is_evaluable(call, sample)) {
            std::vector<Allele> alleles; bool has_ref;
            std::tie(alleles, has_ref) = get_called_alleles(call, sample);
            if (has_ref) alleles.erase(std::cbegin(alleles));
            if (!alleles.empty()) {
                const auto sample_allele_support = compute_allele_support(alleles, assignments, sample);
                for (const auto& p : sample_allele_support) {
                    const auto median_bq = median_base_quality(p.second, p.first);
                    if (median_bq) {
                        if (sample_result) {
                            sample_result = std::min(*sample_result, static_cast<int>(*median_bq));
                        } else {
                            sample_result = *median_bq;
                        }
                    }
                }
            }
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality MedianBaseQuality::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& MedianBaseQuality::do_name() const
{
    return name_;
}

std::string MedianBaseQuality::do_describe() const
{
    return "Median base quality of reads supporting each allele";
}

std::vector<std::string> MedianBaseQuality::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
