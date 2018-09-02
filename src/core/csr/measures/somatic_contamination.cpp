// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_contamination.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/range/combine.hpp>

#include "core/types/allele.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/tools/read_assigner.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/append.hpp"
#include "is_somatic.hpp"
#include "../facets/samples.hpp"
#include "../facets/genotypes.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string SomaticContamination::name_ = "SC";

std::unique_ptr<Measure> SomaticContamination::do_clone() const
{
    return std::make_unique<SomaticContamination>(*this);
}

namespace {

template <typename Container>
void sort_unique(Container& values)
{
    std::sort(std::begin(values), std::end(values));
    values.erase(std::unique(std::begin(values), std::end(values)), std::end(values));
}

auto get_somatic_alleles(const VcfRecord& somatic, const std::vector<SampleName>& somatic_samples,
                         const std::vector<SampleName>& normal_samples)
{
    std::vector<Allele> somatic_sample_alleles {}, normal_sample_alleles {};
    for (const auto& sample : somatic_samples) {
        utils::append(get_called_alleles(somatic, sample, true).first, somatic_sample_alleles);
    }
    for (const auto& sample : normal_samples) {
        utils::append(get_called_alleles(somatic, sample, true).first, normal_sample_alleles);
    }
    sort_unique(somatic_sample_alleles); sort_unique(normal_sample_alleles);
    std::vector<Allele> result {};
    result.reserve(somatic_sample_alleles.size());
    std::set_difference(std::cbegin(somatic_sample_alleles), std::cend(somatic_sample_alleles),
                        std::cbegin(normal_sample_alleles), std::cend(normal_sample_alleles),
                        std::back_inserter(result));
    return result;
}

auto get_somatic_haplotypes(const Facet::GenotypeMap& genotypes, const std::vector<Allele>& somatics)
{
    std::vector<Haplotype> result {};
    if (!somatics.empty()) {
        const auto allele_region = somatics.front().mapped_region();
        for (const auto& p :genotypes) {
            const auto& overlapped_genotypes = overlap_range(p.second, allele_region);
            if (size(overlapped_genotypes) == 1) {
                const auto& genotype = overlapped_genotypes.front();
                for (const auto& haplotype : genotype) {
                    if (std::any_of(std::cbegin(somatics), std::cend(somatics),
                                    [&] (const auto& somatic) { return haplotype.includes(somatic); })) {
                        result.push_back(haplotype);
                    }
                }
            }
        }
        sort_unique(result);
    }
    return result;
}

auto get_somatic_haplotypes(const VcfRecord& somatic, const Facet::GenotypeMap& genotypes,
                            const std::vector<SampleName>& somatic_samples, const std::vector<SampleName>& normal_samples)
{
    const auto somatic_alleles = get_somatic_alleles(somatic, somatic_samples, normal_samples);
    return get_somatic_haplotypes(genotypes, somatic_alleles);
}

template <typename MappableType>
auto copy_overlapped_to_vector(const AmbiguousReadList& reads, const MappableType& mappable)
{
    const auto overlapped = overlap_range(reads, mappable);
    std::vector<AlignedRead> result {};
    result.reserve(size(overlapped));
    std::transform(std::cbegin(overlapped), std::cend(overlapped), std::back_inserter(result),
                   [] (const auto& read) { return read.read; });
    return result;
}

} // namespace

Measure::ResultType SomaticContamination::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    boost::optional<int> result {};
    if (is_somatic(call)) {
        result = 0;
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        const auto somatic_status = boost::get<std::vector<bool>>(IsSomatic(true).evaluate(call, facets));
        std::vector<SampleName> somatic_samples {}, normal_samples {};
        somatic_samples.reserve(samples.size()); normal_samples.reserve(samples.size());
        for (auto tup : boost::combine(samples, somatic_status)) {
            if (tup.get<1>()) {
                somatic_samples.push_back(tup.get<0>());
            } else {
                normal_samples.push_back(tup.get<0>());
            }
        }
        const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
        const auto somatic_haplotypes = get_somatic_haplotypes(call, genotypes, somatic_samples, normal_samples);
        const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
        Genotype<Haplotype> somatic_genotype {static_cast<unsigned>(somatic_haplotypes.size() + 1)};
        HaplotypeProbabilityMap haplotype_priors {};
        haplotype_priors.reserve(somatic_haplotypes.size() + 1);
        for (const auto& haplotype : somatic_haplotypes) {
            somatic_genotype.emplace(haplotype);
            haplotype_priors[haplotype] = -1;
        }
        for (const auto& sample : normal_samples) {
            for (const auto& p : assignments.support.at(sample)) {
                const auto overlapped_reads = copy_overlapped(p.second, call);
                if (!overlapped_reads.empty()) {
                    const Haplotype& assigned_haplotype {p.first};
                    if (!somatic_genotype.contains(assigned_haplotype)) {
                        auto dummy = somatic_genotype;
                        dummy.emplace(assigned_haplotype);
                        haplotype_priors[assigned_haplotype] = 0;
                        const auto support = compute_haplotype_support(dummy, overlapped_reads, haplotype_priors);
                        haplotype_priors.erase(assigned_haplotype);
                        for (const auto& somatic : somatic_haplotypes) {
                            if (support.count(somatic) == 1) {
                                *result += support.at(somatic).size();
                            }
                        }
                    } else {
                        // This could happen if we don't call all 'somatic' alleles on the called somatic haplotype.
                        *result += overlapped_reads.size();
                    }
                }
            }
            if (assignments.ambiguous.count(sample) == 1 && !assignments.ambiguous.at(sample).empty()) {
                const auto ambiguous_reads = copy_overlapped_to_vector(assignments.ambiguous.at(sample), call);
                if (!ambiguous_reads.empty()) {
                    const auto overlapped_genotypes = overlap_range(genotypes.at(sample), call);
                    if (size(overlapped_genotypes) == 1) {
                        const auto& called_genotype = overlapped_genotypes.front();
                        auto dummy = called_genotype;
                        for (const auto& haplotype : somatic_haplotypes) {
                            dummy.emplace(haplotype);
                        }
                        for (const auto& haplptype : called_genotype) {
                            haplotype_priors[haplptype] = 0;
                        }
                        const auto support = compute_haplotype_support(dummy, ambiguous_reads, haplotype_priors);
                        for (const auto& haplotype : called_genotype) {
                            haplotype_priors.erase(haplotype);
                        }
                        for (const auto& somatic : somatic_haplotypes) {
                            if (support.count(somatic) == 1) {
                                *result += support.at(somatic).size();
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality SomaticContamination::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& SomaticContamination::do_name() const
{
    return name_;
}

std::string SomaticContamination::do_describe() const
{
    return "Number of reads supporting a somatic haplotype in the normal";
}

std::vector<std::string> SomaticContamination::do_requirements() const
{
    std::vector<std::string> result {"Samples", "Genotypes", "ReadAssignments"};
    utils::append(IsSomatic(true).requirements(), result);
    sort_unique(result);
    return result;
}

} // namespace csr
} // namespace octopus
