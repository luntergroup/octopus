// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "normal_contamination.hpp"

#include <algorithm>
#include <iterator>
#include <set>

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
#include "../facets/alleles.hpp"
#include "../facets/genotypes.hpp"
#include "../facets/overlapping_reads.hpp"
#include "../facets/pedigree.hpp"

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

} // namespace

Measure::ResultType NormalContamination::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
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
    }
    return result;
}

Measure::ResultCardinality NormalContamination::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
}

const std::string& NormalContamination::do_name() const
{
    return name_;
}

std::string NormalContamination::do_describe() const
{
    return "Fraction of reads that support a called haplotype in another sample with each uncalled allele";
}

std::vector<std::string> NormalContamination::do_requirements() const
{
    std::vector<std::string> result {"Samples", "Alleles", "Genotypes", "OverlappingReads"};
    utils::append(IsSomatic(true).requirements(), result);
    return result;
}

} // namespace csr
} // namespace octopus
