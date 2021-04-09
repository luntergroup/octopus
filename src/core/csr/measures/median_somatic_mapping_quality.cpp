// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "median_somatic_mapping_quality.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/range/combine.hpp>

#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/tools/read_assigner.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/append.hpp"
#include "utils/maths.hpp"
#include "is_somatic.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/genotypes.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MedianSomaticMappingQuality::name_ = "SMQ";

std::unique_ptr<Measure> MedianSomaticMappingQuality::do_clone() const
{
    return std::make_unique<MedianSomaticMappingQuality>(*this);
}

Measure::ValueType MedianSomaticMappingQuality::get_value_type() const
{
    return int {};
}

namespace {

template <typename Container>
void sort_unique(Container& values)
{
    std::sort(std::begin(values), std::end(values));
    values.erase(std::unique(std::begin(values), std::end(values)), std::end(values));
}

auto get_somatic_alleles(const VcfRecord& somatic, 
                         const std::vector<SampleName>& somatic_samples,
                         const std::vector<SampleName>& normal_samples,
                         const Facet::AlleleMap& alleles)
{
    std::vector<Allele> somatic_sample_alleles {}, normal_sample_alleles {};
    for (const auto& sample : somatic_samples) {
        utils::append(get_called(alleles, somatic, sample), somatic_sample_alleles);
    }
    for (const auto& sample : normal_samples) {
        utils::append(get_called(alleles, somatic, sample), normal_sample_alleles);
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
        for (const auto& p : genotypes) {
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

auto get_somatic_haplotypes(const VcfRecord& somatic, 
                            const Facet::GenotypeMap& genotypes,
                            const std::vector<SampleName>& somatic_samples, 
                            const std::vector<SampleName>& normal_samples,
                            const Facet::AlleleMap& alleles)
{
    const auto somatic_alleles = get_somatic_alleles(somatic, somatic_samples, normal_samples, alleles);
    return get_somatic_haplotypes(genotypes, somatic_alleles);
}

} // namespace

Measure::ResultType MedianSomaticMappingQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    Array<Optional<ValueType>> result(samples.size());
    if (is_somatic(call)) {
        const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
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
        if (somatic_samples.empty() || normal_samples.empty()) return result;
        const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
        const auto somatic_haplotypes = get_somatic_haplotypes(call, genotypes, somatic_samples, normal_samples, alleles);
        if (somatic_haplotypes.empty()) return result;
        const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).haplotypes;
        for (std::size_t s {0}; s < samples.size(); ++s) {
            const auto& sample = samples[s];
            if (normal_samples.empty()
                || std::find(std::cbegin(somatic_samples), std::cend(somatic_samples), sample) != std::cend(somatic_samples)) {
                std::vector<AlignedRead::MappingQuality> somatic_mqs {};
                for (const auto& haplotype : somatic_haplotypes) {
                    if (assignments.at(sample).assigned_wrt_reference.count(haplotype) == 0) {
                        // TODO: This should never happen...
                        somatic_mqs.clear();
                        break;
                    }
                    const auto& somatic_support = assignments.at(sample).assigned_wrt_reference.at(haplotype);
                    somatic_mqs.reserve(somatic_mqs.size() + somatic_support.size());
                    std::transform(std::cbegin(somatic_support), std::cend(somatic_support), std::back_inserter(somatic_mqs),
                                   [] (const AlignedRead& read) { return read.mapping_quality(); });
                }
                if (!somatic_mqs.empty()) {
                    result[s] = static_cast<int>(maths::median(somatic_mqs));
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality MedianSomaticMappingQuality::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
}

const std::string& MedianSomaticMappingQuality::do_name() const
{
    return name_;
}

std::string MedianSomaticMappingQuality::do_describe() const
{
    return "Median mapping quality of reads assigned to called somatic haplotypes";
}

std::vector<std::string> MedianSomaticMappingQuality::do_requirements() const
{
    std::vector<std::string> result {"Samples", "Alleles", "Genotypes", "ReadAssignments"};
    utils::append(IsSomatic(true).requirements(), result);
    sort_unique(result);
    return result;
}
    
} // namespace csr
} // namespace octopus
