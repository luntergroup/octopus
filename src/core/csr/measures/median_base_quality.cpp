// Copyright (c) 2015-2020 Daniel Cooke
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
#include "utils/maths.hpp"
#include "utils/append.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MedianBaseQuality::name_ = "BQ";

std::unique_ptr<Measure> MedianBaseQuality::do_clone() const
{
    return std::make_unique<MedianBaseQuality>(*this);
}

Measure::ValueType MedianBaseQuality::get_value_type() const
{
    return int {};
}

namespace {

bool is_canonical(const VcfRecord::NucleotideSequence& allele) noexcept
{
    return !(allele == vcfspec::missingValue || allele == vcfspec::deleteMaskAllele);
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
    Measure::Optional<Measure::ValueType> result {};
    if (!is_indel(allele)) {
        std::vector<AlignedRead::BaseQuality> base_qualities {};
        base_qualities.reserve(reads.size() * sequence_size(allele));
        for (const auto& read : reads) {
            if (overlaps(read.get(), allele)) {
                utils::append(copy_base_qualities(read, mapped_region(allele)), base_qualities);
            }
        }
        if (!base_qualities.empty()) {
            result = static_cast<int>(maths::median(base_qualities));
        }
    }
    return result;
}

} // namespace

Measure::ResultType MedianBaseQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    const auto num_alleles = call.num_alt() + 1;
    Array<Array<Optional<ValueType>>> result(samples.size(), Array<Optional<ValueType>>(num_alleles));
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample_alleles = get(alleles, call, samples[s]);
        assert(sample_alleles.size() == num_alleles);
        const auto& support = assignments.at(samples[s]);
        for (std::size_t a {0}; a < num_alleles; ++a) {
            if (sample_alleles[a]) {
                const auto& allele = *sample_alleles[a];
                const auto support_set_itr = support.find(allele);
                if (support_set_itr != std::cend(support)) {
                    result[s][a] = median_base_quality(support_set_itr->second, allele);
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality MedianBaseQuality::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
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
    return {"Samples", "Alleles", "ReadAssignments"};
}

boost::optional<Measure::Aggregator> MedianBaseQuality::do_aggregator() const noexcept
{
    return Measure::Aggregator::min_tail;
}

} // namespace csr
} // namespace octopus
