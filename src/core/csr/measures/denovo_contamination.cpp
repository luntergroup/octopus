// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_contamination.hpp"

#include <algorithm>
#include <iterator>
#include <array>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/range/combine.hpp>

#include "basics/trio.hpp"
#include "core/types/allele.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/tools/read_assigner.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/genotype_reader.hpp"
#include "utils/concat.hpp"
#include "utils/append.hpp"
#include "is_denovo.hpp"
#include "../facets/samples.hpp"
#include "../facets/genotypes.hpp"
#include "../facets/pedigree.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string DeNovoContamination::name_ = "DC";

std::unique_ptr<Measure> DeNovoContamination::do_clone() const
{
    return std::make_unique<DeNovoContamination>(*this);
}

namespace {

bool is_denovo(const VcfRecord& call, const Measure::FacetMap& facets)
{
    return boost::get<bool>(IsDenovo(false).evaluate(call, facets));
}

auto find_child_idx(const std::vector<SampleName>& samples, const octopus::Pedigree& pedigree)
{
    assert(samples.size() == 3);
    if (is_parent_of(samples[0], samples[1], pedigree)) {
        return 1;
    } else if (is_parent_of(samples[1], samples[0], pedigree)) {
        return 0;
    } else {
        return 2;
    }
}

template <typename Container>
void sort_unique(Container& values)
{
    std::sort(std::begin(values), std::end(values));
    values.erase(std::unique(std::begin(values), std::end(values)), std::end(values));
}

auto get_denovo_alleles(const VcfRecord& denovo, const Trio& trio)
{
    auto parent_alleles = concat(get_called_alleles(denovo, trio.mother(), true).first,
                                 get_called_alleles(denovo, trio.father(), true).first);
    auto child_alleles = get_called_alleles(denovo, trio.child(), true).first;
    sort_unique(parent_alleles); sort_unique(child_alleles);
    std::vector<Allele> result {};
    result.reserve(child_alleles.size());
    std::set_difference(std::cbegin(child_alleles), std::cend(child_alleles),
                        std::cbegin(parent_alleles), std::cend(parent_alleles),
                        std::back_inserter(result));
    return result;
}

auto get_denovo_haplotypes(const Facet::GenotypeMap& genotypes, const std::vector<Allele>& denovos)
{
    std::vector<Haplotype> result {};
    if (!denovos.empty()) {
        const auto allele_region = denovos.front().mapped_region();
        for (const auto& p :genotypes) {
            const auto& overlapped_genotypes = overlap_range(p.second, allele_region);
            if (size(overlapped_genotypes) == 1) {
                const auto& genotype = overlapped_genotypes.front();
                for (const auto& haplotype : genotype) {
                    if (std::any_of(std::cbegin(denovos), std::cend(denovos),
                                    [&] (const auto& denovo) { return haplotype.includes(denovo); })) {
                        result.push_back(haplotype);
                    }
                }
            }
        }
        sort_unique(result);
    }
    return result;
}

auto get_denovo_haplotypes(const VcfRecord& denovo, const Facet::GenotypeMap& genotypes, const Trio& trio)
{
    const auto denovo_alleles = get_denovo_alleles(denovo, trio);
    return get_denovo_haplotypes(genotypes, denovo_alleles);
}

template <typename Container, typename MappableType>
void remove_not_contained(Container& reads, const MappableType& mappable)
{
    reads.erase(std::remove_if(std::begin(reads), std::end(reads),
                               [&mappable] (const auto& read) { return !contains(mappable, read); }),
                std::end(reads));
}

} // namespace

Measure::ResultType DeNovoContamination::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    boost::optional<int> result {};
    if (is_denovo(call, facets)) {
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        const auto& pedigree = get_value<Pedigree>(facets.at("Pedigree"));
        assert(is_trio(samples, pedigree)); // TODO: Implement for general pedigree
        const auto trio = *make_trio(samples[find_child_idx(samples, pedigree)], pedigree);
        const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
        const auto denovo_haplotypes = get_denovo_haplotypes(call, genotypes, trio);
        if (denovo_haplotypes.empty()) {
            return result;
        }
        const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).support;
        Genotype<Haplotype> denovo_genotype {static_cast<unsigned>(denovo_haplotypes.size() + 1)};
        HaplotypeProbabilityMap haplotype_priors {};
        haplotype_priors.reserve(denovo_haplotypes.size() + 1);
        for (const auto& haplotype : denovo_haplotypes) {
            denovo_genotype.emplace(haplotype);
            haplotype_priors[haplotype] = -1;
        }
        const std::array<SampleName, 2> parents {trio.mother(), trio.father()};
        result = 0;
        for (const auto& sample : parents) {
            for (const auto& p : assignments.at(sample)) {
                auto supporting_reads = copy_overlapped(p.second, call);
                if (!supporting_reads.empty()) {
                    const Haplotype& assigned_haplotype {p.first};
                    if (!denovo_genotype.contains(assigned_haplotype)) {
                        if (!is_same_region(assigned_haplotype, denovo_genotype)) {
                            remove_not_contained(supporting_reads, assigned_haplotype);
                        }
                        if (!supporting_reads.empty()) {
                            auto dummy = denovo_genotype;
                            dummy.emplace(assigned_haplotype);
                            haplotype_priors[assigned_haplotype] = 0;
                            const auto support = compute_haplotype_support(dummy, supporting_reads, haplotype_priors);
                            haplotype_priors.erase(assigned_haplotype);
                            for (const auto& denovo : denovo_haplotypes) {
                                if (support.count(denovo) == 1) {
                                    *result += support.at(denovo).size();
                                }
                            }
                        }
                    } else {
                        // This could happen if we don't call all 'de novo' alleles on the called de novo haplotype.
                        *result += supporting_reads.size();
                    }
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality DeNovoContamination::do_cardinality() const noexcept
{
    return ResultCardinality::one;
}

const std::string& DeNovoContamination::do_name() const
{
    return name_;
}

std::string DeNovoContamination::do_describe() const
{
    return "Number of reads supporting a de novo haplotype in the normal";
}

std::vector<std::string> DeNovoContamination::do_requirements() const
{
    std::vector<std::string> result {"Samples", "Genotypes", "ReadAssignments", "Pedigree"};
    utils::append(IsDenovo(false).requirements(), result);
    sort_unique(result);
    return result;
}

} // namespace csr
} // namespace octopus
