// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_assignments.hpp"

#include "core/tools/read_realigner.hpp"
#include "utils/genotype_reader.hpp"

namespace octopus { namespace csr {

const std::string ReadAssignments::name_ {"ReadAssignments"};

namespace {

template <typename Mappable>
auto copy_overlapped_to_vector(const ReadContainer& reads, const Mappable& mappable)
{
    const auto overlapped = overlap_range(reads, mappable);
    return std::vector<AlignedRead> {std::cbegin(overlapped), std::cend(overlapped)};
}

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles, 
                       const Facet::SampleSupportMap& support, 
                       const Facet::SampleAmbiguityMap& ambiguous,
                       const SampleName& sample)
{
    if (support.count(sample) == 1) {
        if (ambiguous.count(sample) == 1) {
            return compute_allele_support(alleles, support.at(sample), ambiguous.at(sample));
        } else {
            return compute_allele_support(alleles, support.at(sample));
        }
    } else {
        if (ambiguous.count(sample) == 1) {
            const HaplotypeSupportMap empty {};
            return compute_allele_support(alleles, empty, ambiguous.at(sample));
        } else {
            return {};
        }
    }
}

} // namespace

ReadAssignments::ReadAssignments(const ReferenceGenome& reference,
                                 const GenotypeMap& genotypes,
                                 const ReadMap& reads,
                                 const std::vector<VcfRecord>& calls)
: ReadAssignments {reference, genotypes, reads, calls, {}} {}

ReadAssignments::ReadAssignments(const ReferenceGenome& reference,
                                 const GenotypeMap& genotypes,
                                 const ReadMap& reads,
                                 const std::vector<VcfRecord>& calls,
                                 HaplotypeLikelihoodModel model)
: result_ {}
, likelihood_model_ {std::move(model)}
{
    const auto num_samples = genotypes.size();
    result_.support.reserve(num_samples);
    result_.ambiguous.reserve(num_samples);
    for (const auto& p : genotypes) {
        const auto& sample = p.first;
        const auto& sample_genotypes = p.second;
        result_.support[sample].reserve(sample_genotypes.size());
        for (const auto& genotype : sample_genotypes) {
            auto local_reads = copy_overlapped_to_vector(reads.at(sample), genotype);
            for (const auto& haplotype : genotype) {
                // So every called haplotype appears in support map, even if no read support
                result_.support[sample][haplotype] = {};
            }
            if (!local_reads.empty()) {
                HaplotypeSupportMap genotype_support {};
                if (!is_homozygous(genotype)) {
                    genotype_support = compute_haplotype_support(genotype, local_reads, result_.ambiguous[sample], likelihood_model_);
                } else {
                    if (is_reference(genotype[0])) {
                        genotype_support[genotype[0]] = std::move(local_reads);
                    } else {
                        auto augmented_genotype = genotype;
                        Haplotype ref {mapped_region(genotype), reference};
                        result_.support[sample][ref] = {};
                        augmented_genotype.emplace(std::move(ref));
                        genotype_support = compute_haplotype_support(augmented_genotype, local_reads, result_.ambiguous[sample], likelihood_model_);
                    }
                }
                for (auto& s : genotype_support) {
                    safe_realign_to_reference(s.second, s.first, likelihood_model_);
                    std::sort(std::begin(s.second), std::end(s.second));
                    result_.support[sample][s.first] = std::move(s.second);
                }
            }
        }
        for (const auto& call : calls) {
            auto alleles = get_called_alleles(call, sample).first;
            auto allele_support = compute_allele_support(alleles, result_.support, result_.ambiguous, sample);
            for (auto& allele : alleles) {
                result_.alleles[sample].emplace(std::move(allele), std::move(allele_support.at(allele)));
            }
        }
    }
}

Facet::ResultType ReadAssignments::do_get() const
{
    return std::cref(result_);
}

AlleleSupportMap
compute_allele_support(const std::vector<Allele>& alleles, const Facet::SupportMaps& assignments, const SampleName& sample)
{
    if (assignments.support.count(sample) == 1) {
        if (assignments.ambiguous.count(sample) == 1) {
            return compute_allele_support(alleles, assignments.support.at(sample), assignments.ambiguous.at(sample));
        } else {
            return compute_allele_support(alleles, assignments.support.at(sample));
        }
    } else {
        if (assignments.ambiguous.count(sample) == 1) {
            const HaplotypeSupportMap empty {};
            return compute_allele_support(alleles, empty, assignments.ambiguous.at(sample));
        } else {
            return {};
        }
    }
}

} // namespace csr
} // namespace octopus
