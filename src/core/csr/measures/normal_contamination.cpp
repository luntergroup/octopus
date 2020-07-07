// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "normal_contamination.hpp"

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

const std::string NormalContamination::name_ = "NC";

std::unique_ptr<Measure> NormalContamination::do_clone() const
{
    return std::make_unique<NormalContamination>(*this);
}

Measure::ValueType NormalContamination::get_value_type() const
{
    return double {};
}

namespace {

auto extract_called_alleles(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    return get_called_alleles(call, sample, ReferencePadPolicy::trim_alt_alleles).first;
}

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
        utils::append(extract_called_alleles(somatic, sample), somatic_sample_alleles);
    }
    for (const auto& sample : normal_samples) {
        utils::append(extract_called_alleles(somatic, sample), normal_sample_alleles);
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

Measure::ResultType NormalContamination::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    Optional<ValueType> result {};
    if (is_somatic(call)) {
        std::size_t contamination {0}, total_overlapped {0};
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        const auto somatic_status = boost::get<Array<ValueType>>(IsSomatic(true).evaluate(call, facets));
        std::vector<SampleName> somatic_samples {}, normal_samples {};
        somatic_samples.reserve(samples.size()); normal_samples.reserve(samples.size());
        for (auto tup : boost::combine(samples, somatic_status)) {
            if (boost::get<bool>(tup.get<1>())) {
                somatic_samples.push_back(tup.get<0>());
            } else {
                normal_samples.push_back(tup.get<0>());
            }
        }
        if (normal_samples.empty()) return result;
        const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
        const auto somatic_haplotypes = get_somatic_haplotypes(call, genotypes, somatic_samples, normal_samples);
        const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).haplotypes;
        Genotype<Haplotype> somatic_genotype {static_cast<unsigned>(somatic_haplotypes.size() + 1)};
        for (const auto& haplotype : somatic_haplotypes) {
            somatic_genotype.emplace(haplotype);
        }
        for (const auto& sample : normal_samples) {
            HaplotypeProbabilityMap haplotype_priors {};
            haplotype_priors.reserve(somatic_haplotypes.size() + 1);
            for (const auto& haplotype : somatic_haplotypes) {
                haplotype_priors[haplotype] = -1;
            }
            for (const auto& p : assignments.at(sample).assigned_wrt_reference) {
                const auto overlapped_reads = copy_overlapped(p.second, call);
                if (!overlapped_reads.empty()) {
                    const Haplotype& assigned_haplotype {p.first};
                    if (!contains(somatic_genotype, assigned_haplotype)) {
                        auto dummy = somatic_genotype;
                        dummy.emplace(assigned_haplotype);
                        haplotype_priors[assigned_haplotype] = 0;
                        const auto support = compute_haplotype_support(dummy, overlapped_reads, haplotype_priors);
                        haplotype_priors.erase(assigned_haplotype);
                        for (const auto& somatic : somatic_haplotypes) {
                            if (support.count(somatic) == 1) {
                                contamination += support.at(somatic).size();
                            }
                        }
                    } else {
                        // This could happen if we don't call all 'somatic' alleles on the called somatic haplotype.
                        contamination += overlapped_reads.size();
                    }
                    total_overlapped += overlapped_reads.size();
                }
            }
            if (!assignments.at(sample).ambiguous_wrt_reference.empty()) {
                const auto ambiguous_reads = copy_overlapped_to_vector(assignments.at(sample).ambiguous_wrt_reference, call);
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
                        for (const auto& somatic : somatic_haplotypes) {
                            if (support.count(somatic) == 1) {
                                contamination += support.at(somatic).size();
                            }
                        }
                    }
                    total_overlapped += ambiguous_reads.size();
                }
            }
        }
        result = total_overlapped > 0 ? static_cast<double>(contamination) / total_overlapped : 0.0;
    }
    return result;
}

Measure::ResultCardinality NormalContamination::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& NormalContamination::do_name() const
{
    return name_;
}

std::string NormalContamination::do_describe() const
{
    return "Fraction of overlapping reads supporting a somatic haplotype in the normal";
}

std::vector<std::string> NormalContamination::do_requirements() const
{
    std::vector<std::string> result {"Samples", "Genotypes", "ReadAssignments"};
    utils::append(IsSomatic(true).requirements(), result);
    sort_unique(result);
    return result;
}

} // namespace csr
} // namespace octopus
