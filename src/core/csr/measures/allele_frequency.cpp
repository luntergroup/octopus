// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_frequency.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string AlleleFrequency::name_ = "AF";

std::unique_ptr<Measure> AlleleFrequency::do_clone() const
{
    return std::make_unique<AlleleFrequency>(*this);
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

} // namespace

Measure::ResultType AlleleFrequency::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).support;
    std::vector<boost::optional<double>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        boost::optional<double> sample_result {};
        if (is_evaluable(call, sample)) {
            const auto& sample_assignments = assignments.at(sample);
            std::vector<Allele> alleles; bool has_ref;
            std::tie(alleles, has_ref) = get_called_alleles(call, sample, true);
            assert(!alleles.empty());
            std::size_t read_count {0};
            std::vector<unsigned> allele_counts(alleles.size());
            for (const auto& p : sample_assignments) {
                const auto& haplotype = p.first;
                const auto& reads = p.second;
                const auto haplotype_support_depth = count_overlapped(reads, call);
                if (haplotype_support_depth > 0) {
                    std::transform(std::cbegin(alleles), std::cend(alleles), std::cbegin(allele_counts), std::begin(allele_counts),
                                   [&] (const auto& allele, auto count) {
                                       if (haplotype.includes(allele)) {
                                           count += haplotype_support_depth;
                                       }
                                       return count;
                                   });
                    read_count += haplotype_support_depth;
                }
            }
            if (read_count > 0) {
                auto first_called_count_itr = std::cbegin(allele_counts);
                if (has_ref) ++first_called_count_itr;
                assert(first_called_count_itr != std::cend(allele_counts));
                const auto min_count_itr = std::min_element(first_called_count_itr, std::cend(allele_counts));
                sample_result = static_cast<double>(*min_count_itr) / read_count;
            }
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality AlleleFrequency::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& AlleleFrequency::do_name() const
{
    return name_;
}

std::string AlleleFrequency::do_describe() const
{
    return "Minor allele frequency of ALT alleles";
}

std::vector<std::string> AlleleFrequency::do_requirements() const
{
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
