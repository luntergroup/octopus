// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "contamination.hpp"

#include <algorithm>
#include <iterator>
#include <set>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "core/types/allele.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/tools/read_assigner.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/append.hpp"
#include "is_somatic.hpp"
#include "is_denovo.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/genotypes.hpp"
#include "../facets/overlapping_reads.hpp"
#include "../facets/pedigree.hpp"

namespace octopus { namespace csr {

const std::string Contamination::name_ = "CF";

std::unique_ptr<Measure> Contamination::do_clone() const
{
    return std::make_unique<Contamination>(*this);
}

Measure::ValueType Contamination::get_value_type() const
{
    return double {};
}

namespace {

bool is_denovo(const VcfRecord& call, const Measure::FacetMap& facets)
{
    return get_value_type<bool>(IsDenovo(false).evaluate(call, facets));
}

auto convert_to_bools(const Measure::Array<Measure::ValueType>& values)
{
    std::vector<bool> result(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::begin(result),
                   [] (auto value) { return boost::get<bool>(value); });
    return result;
}

const auto& get_genotype(const Facet::GenotypeMap& genotypes, const VcfRecord& call, const SampleName& sample)
{
    return overlap_range(genotypes.at(sample), call).front();
}

auto convert_to_bools(const std::vector<std::string>& values)
{
    std::vector<bool> result(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::begin(result),
                   [] (const auto& value) { return value == "1" ? true : false; });
    return result;
}

auto get_haplotype_somatic_statuses(const VcfRecord& call, const SampleName& sample)
{
    return convert_to_bools(call.get_sample_value(sample, "HSS"));
}

template <typename Range, typename MappableType>
auto copy_overlapped_to_vector(const Range& mappables, const MappableType& mappable)
{
    const auto overlapped = overlap_range(mappables, mappable);
    return std::vector<typename Range::value_type> {std::cbegin(overlapped), std::cend(overlapped)};
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

} // namespace

Measure::ResultType Contamination::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto num_alleles = 1 + call.num_alt();
    Array<Array<Optional<ValueType>>> result(samples.size(), Array<Optional<ValueType>>(num_alleles));
    if (is_somatic(call)) {
        const auto somatic_status = convert_to_bools(boost::get<Array<ValueType>>(IsSomatic(true).evaluate(call, facets)));
        std::vector<std::size_t> somatic_sample_indices {}, non_somatic_sample_indices {};
        for (std::size_t s {0}; s < samples.size(); ++s) {
            if (somatic_status[s]) {
                somatic_sample_indices.push_back(s);
            } else {
                non_somatic_sample_indices.push_back(s);
            }
        }
        if (!non_somatic_sample_indices.empty()) {
            const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
            std::set<Haplotype> somatic_haplotypes {};
            for (const auto sample_idx : somatic_sample_indices) {
                const auto& sample = samples[sample_idx];
                const auto haplotype_somatic_status = get_haplotype_somatic_statuses(call, sample);
                const auto& genotype = get_genotype(genotypes, call, sample);
                assert(genotype.ploidy() == call.ploidy(sample));
                for (unsigned p {0}; p < genotype.ploidy(); ++p) {
                    if (haplotype_somatic_status[p]) {
                        somatic_haplotypes.insert(genotype[p]);
                    }
                }
            }
            assert(somatic_haplotypes.size() > 0);
            // All non-somatic samples are assumed to have the same genotype - the germline
            const auto& germline_genotype = get_genotype(genotypes, call, samples[non_somatic_sample_indices.front()]);
            Genotype<Haplotype> contaminated_germline_genotype {germline_genotype.ploidy() + static_cast<unsigned>(somatic_haplotypes.size())};
            for (const auto& haplotype : germline_genotype) {
                contaminated_germline_genotype.emplace(haplotype);
            }
            for (const auto& haplotype : somatic_haplotypes) {
                contaminated_germline_genotype.emplace(haplotype);
            }
            // Put weak priors against the somatic haplotypes to favour germline for ambiguous reads
            HaplotypeProbabilityMap haplotype_priors {};
            haplotype_priors.reserve(contaminated_germline_genotype.ploidy());
            for (const auto& haplotype : germline_genotype) {
                haplotype_priors[haplotype] = 0;
            }
            for (const auto& haplotype : somatic_haplotypes) {
                haplotype_priors[haplotype] = -0.5;
            }
            const auto& reads = get_value<OverlappingReads>(facets.at("OverlappingReads"));
            for (const auto sample_idx : non_somatic_sample_indices) {
                const auto& sample = samples[sample_idx];
                const auto sample_reads = copy_overlapped_to_vector(reads.at(sample), call);
                if (!sample_reads.empty()) {
                    const auto support = compute_haplotype_support(contaminated_germline_genotype, sample_reads, haplotype_priors);
                    const auto& sample_alleles = get(alleles, call, sample);
                    assert(sample_alleles.size() == num_alleles);
                    for (std::size_t a {0}; a < num_alleles; ++a) {
                        if (!sample_alleles[a]) {
                            const auto is_allele_called = [&] (auto sample_idx) { return get(alleles, call, samples[sample_idx])[a].has_value(); };
                            const auto somatic_sample_itr = std::find_if(std::cbegin(somatic_sample_indices), std::cend(somatic_sample_indices), is_allele_called);
                            if (somatic_sample_itr != std::cend(somatic_sample_indices)) {
                                const auto& somatic_allele = *get(alleles, call, samples[*somatic_sample_itr])[a];
                                unsigned contamination {0};
                                for (const auto& haplotype : somatic_haplotypes) {
                                    if (support.count(haplotype) == 1 && haplotype.contains(somatic_allele)) {
                                        contamination += support.at(haplotype).size();
                                    }
                                }
                                result[sample_idx][a] = static_cast<double>(contamination) / sample_reads.size();
                            }
                        }
                    }
                }
            }
        }
    } else if (is_denovo(call, facets)) {
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

Measure::ResultCardinality Contamination::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
}

const std::string& Contamination::do_name() const
{
    return name_;
}

std::string Contamination::do_describe() const
{
    return "Fraction of reads that support a haplotype called in another sample containing each the uncalled allele";
}

std::vector<std::string> Contamination::do_requirements() const
{
    std::vector<std::string> result {"Samples", "Genotypes", "Alleles", "OverlappingReads", "Pedigree"};
    utils::append(IsSomatic(true).requirements(), result);
    return result;
}

boost::optional<Measure::Aggregator> Contamination::do_aggregator() const noexcept
{
    return Measure::Aggregator::sum;
}

} // namespace csr
} // namespace octopus
