// Copyright (c) 2015-2020 Daniel Cooke
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
#include "utils/concat.hpp"
#include "utils/append.hpp"
#include "is_denovo.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/genotypes.hpp"
#include "../facets/overlapping_reads.hpp"
#include "../facets/pedigree.hpp"

namespace octopus { namespace csr {

const std::string DeNovoContamination::name_ = "DC";

std::unique_ptr<Measure> DeNovoContamination::do_clone() const
{
    return std::make_unique<DeNovoContamination>(*this);
}

Measure::ValueType DeNovoContamination::get_value_type() const
{
    return double {};
}

namespace {

bool is_denovo(const VcfRecord& call, const Measure::FacetMap& facets)
{
    return get_value_type<bool>(IsDenovo(false).evaluate(call, facets));
}

const auto& get_genotype(const Facet::GenotypeMap& genotypes, const VcfRecord& call, const SampleName& sample)
{
    return overlap_range(genotypes.at(sample), call).front();
}

template <typename T>
std::size_t find_index(const T& value, const std::vector<T>& values)
{
    return std::distance(std::cbegin(values), std::find(std::cbegin(values), std::cend(values), value));
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

template <typename Range, typename MappableType>
auto copy_overlapped_to_vector(const Range& mappables, const MappableType& mappable)
{
    const auto overlapped = overlap_range(mappables, mappable);
    return std::vector<typename Range::value_type> {std::cbegin(overlapped), std::cend(overlapped)};
}

} // namespace

Measure::ResultType DeNovoContamination::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto num_alleles = 1 + call.num_alt();
    Array<Array<Optional<ValueType>>> result(samples.size(), Array<Optional<ValueType>>(num_alleles));
    if (is_denovo(call, facets)) {
        const auto& pedigree = get_value<Pedigree>(facets.at("Pedigree"));
        assert(is_trio(samples, pedigree)); // TODO: Implement for general pedigree
        const auto trio = *make_trio(samples[find_child_idx(samples, pedigree)], pedigree);
        const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
        const auto& child_genotype = get_genotype(genotypes, call, trio.child());
        const auto& mother_genotype = get_genotype(genotypes, call, trio.mother());
        const auto& father_genotype = get_genotype(genotypes, call, trio.father());
        std::set<Haplotype> denovo_haplotypes {};
        for (const auto& haplotype : child_genotype) {
            if (!contains(mother_genotype, haplotype) && !contains(father_genotype, haplotype)) {
                denovo_haplotypes.insert(haplotype);
            }
        }
        assert(!denovo_haplotypes.empty());
        const std::array<std::size_t, 2> parent_indices {find_index(trio.mother(), samples), find_index(trio.father(), samples)};
        const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
        HaplotypeProbabilityMap haplotype_priors {};
        // Put weak priors against the denovo haplotypes
        for (const auto& haplotype : denovo_haplotypes) {
            haplotype_priors[haplotype] = -0.5;
        }
        const auto& child_alleles = get(alleles, call, trio.child());
        for (const auto parent_idx : parent_indices) {
            const auto& parent = samples[parent_idx];
            const auto& genotype = parent_idx == parent_indices.front() ? mother_genotype : father_genotype;
            Genotype<Haplotype> contaminated_genotype {genotype.ploidy() + static_cast<unsigned>(denovo_haplotypes.size())};
            for (const auto& haplotype : genotype) {
                contaminated_genotype.emplace(haplotype);
            }
            for (const auto& haplotype : denovo_haplotypes) {
                contaminated_genotype.emplace(haplotype);
            }
            for (const auto& haplotype : genotype) {
                haplotype_priors[haplotype] = 0;
            }
            const auto parent_reads = copy_overlapped_to_vector(reads.at(parent), call);
            if (!parent_reads.empty()) {
                    const auto support = compute_haplotype_support(contaminated_genotype, parent_reads, haplotype_priors);
                    const auto& parent_alleles = get(alleles, call, parent);
                    assert(parent_alleles.size() == num_alleles);
                    for (std::size_t a {0}; a < num_alleles; ++a) {
                        if (!parent_alleles[a] && child_alleles[a]) {
                            unsigned contamination {0};
                            for (const auto& haplotype : denovo_haplotypes) {
                                if (support.count(haplotype) == 1 && haplotype.contains(*child_alleles[a])) {
                                    contamination += support.at(haplotype).size();
                                }
                            }
                            result[parent_idx][a] = static_cast<double>(contamination) / parent_reads.size();
                        }
                    }
            }
        }
    }
    return result;
}

Measure::ResultCardinality DeNovoContamination::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
}

const std::string& DeNovoContamination::do_name() const
{
    return name_;
}

std::string DeNovoContamination::do_describe() const
{
    return "Fraction of reads that support a called haplotype in the child containing each uncalled allele";
}

std::vector<std::string> DeNovoContamination::do_requirements() const
{
    std::vector<std::string> result {"Samples", "Alleles", "Genotypes", "Pedigree"};
    utils::append(IsDenovo(false).requirements(), result);
    return result;
}

} // namespace csr
} // namespace octopus
