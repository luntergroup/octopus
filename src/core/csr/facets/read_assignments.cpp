// Copyright (c) 2015-2021 Daniel Cooke
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
                       const Facet::SupportMaps::HaplotypeSupportMaps& support, 
                       const SampleName& sample,
                       const ReferenceGenome& reference)
{
    const auto overlap_aware_includes = [&] (const Haplotype& haplotype, const Allele& allele) -> bool {
        if (!contains(haplotype, mapped_region(allele))) {
            return false;
        } else if (haplotype.includes(allele)) {
            if (is_position(allele) && is_reference(allele, reference)) {
                // not an insertion padding base
                return haplotype.sequence_size(tail_region(allele)) == 0;
            } else {
                return true;
            }
        } else if (is_empty_region(allele) && is_sequence_empty(allele)) {
            // Spanning deletion
            return false;
        } else if (is_reference(allele, reference)) {
            if (!haplotype.contains(allele)) return false;
            const Allele spanning_deletion {head_region(allele), ""};
            if (haplotype.includes(spanning_deletion)) return false;
            if (begins_equal(allele, haplotype) && is_position(allele)) {
                // not an insertion padding base
                return haplotype.sequence_size(tail_region(allele)) == 0;
            } else {
                return true;
            }
        } else {
            return false;
        }
    };
    return compute_allele_support(alleles, support.assigned_wrt_reference, support.ambiguous_wrt_reference, overlap_aware_includes);
}

template <typename T1, typename T2,
          typename BinaryPredicate = std::less<std::pair<T1, T2>>>
void sort_together(std::vector<T1>& first, std::vector<T2>& second,
                   BinaryPredicate pred = std::less<std::pair<T1, T2>> {})
{
    assert(first.size() == second.size());
    std::vector<std::pair<T1, T2>> zipped(first.size());
    for (std::size_t i {0}; i < first.size(); ++i) {
        zipped[i] = std::make_pair(std::move(first[i]), std::move(second[i]));
    }
    std::sort(std::begin(zipped), std::end(zipped), pred);
    for (std::size_t i {0}; i < first.size(); ++i) {
        first[i]  = std::move(zipped[i].first);
        second[i] = std::move(zipped[i].second);
    }
}

} // namespace

ReadAssignments::ReadAssignments(const ReferenceGenome& reference,
                                 const GenotypeMap& genotypes,
                                 const ReadMap& reads,
                                 const std::vector<VcfRecord>& calls,
                                 OptionalThreadPool workers)
: ReadAssignments {reference, genotypes, reads, calls, {}, workers} {}

ReadAssignments::ReadAssignments(const ReferenceGenome& reference,
                                 const GenotypeMap& genotypes,
                                 const ReadMap& reads,
                                 const std::vector<VcfRecord>& calls,
                                 HaplotypeLikelihoodModel model,
                                 OptionalThreadPool workers)
: result_ {}
, likelihood_model_ {std::move(model)}
{
    const auto num_samples = genotypes.size();
    result_.haplotypes.reserve(num_samples);
    for (const auto& p : genotypes) {
        const auto& sample = p.first;
        const auto& sample_genotypes = p.second;
        result_.haplotypes[sample].assigned_wrt_reference.reserve(sample_genotypes.size());
        result_.haplotypes[sample].assigned_wrt_haplotype.reserve(sample_genotypes.size());
        result_.alleles[sample] = {}; // make sure sample is present
        for (const auto& genotype : sample_genotypes) {
            auto local_reads = copy_overlapped_to_vector(reads.at(sample), genotype);
            for (const auto& haplotype : genotype) {
                // So every called haplotype appears in support map, even if no read support
                result_.haplotypes[sample].assigned_wrt_reference[haplotype] = {};
                result_.haplotypes[sample].assigned_wrt_haplotype[haplotype] = {};
                result_.haplotypes[sample].assigned_likelihoods[haplotype] = {};
            }
            if (!local_reads.empty()) {
                const AssignmentConfig config {};
                // Try to assign each read to a haplotype
                HaplotypeSupportMap genotype_support {};
                if (is_heterozygous(genotype)) {
                    genotype_support = compute_haplotype_support(genotype, local_reads, result_.haplotypes[sample].ambiguous_wrt_haplotype, likelihood_model_, config, workers);
                } else {
                    if (is_reference(genotype[0])) {
                        genotype_support[genotype[0]] = std::move(local_reads);
                    } else {
                        auto augmented_genotype = genotype;
                        Haplotype ref {mapped_region(genotype), reference};
                        result_.haplotypes[sample].assigned_wrt_reference[ref] = {};
                        augmented_genotype.emplace(std::move(ref));
                        genotype_support = compute_haplotype_support(augmented_genotype, local_reads, result_.haplotypes[sample].ambiguous_wrt_haplotype, likelihood_model_, config, workers);
                    }
                }
                // Realign assigned reads
                for (auto& s : genotype_support) {
                    const Haplotype& haplotype {s.first};
                    auto& assigned_reads = s.second;
                    auto& likelihoods = result_.haplotypes[sample].assigned_likelihoods[haplotype];
                    safe_realign(assigned_reads, haplotype, likelihood_model_, likelihoods, workers);
                    sort_together(assigned_reads, likelihoods);
                    result_.haplotypes[sample].assigned_wrt_haplotype[haplotype] = assigned_reads;
                    rebase(assigned_reads, haplotype, workers);
                    std::sort(std::begin(assigned_reads), std::end(assigned_reads));
                    result_.haplotypes[sample].assigned_wrt_reference[haplotype] = std::move(assigned_reads);
                }
                // Realign ambiguous reads
                std::unordered_map<Haplotype, std::vector<std::size_t>> possible_ambiguous_assignments {};
                auto& ambiguous_reads = result_.haplotypes[sample].ambiguous_wrt_haplotype;
                for (std::size_t ambiguous_read_idx {0}; ambiguous_read_idx < ambiguous_reads.size(); ++ambiguous_read_idx) {
                    const auto& read = ambiguous_reads[ambiguous_read_idx];
                    if (read.haplotypes) {
                        possible_ambiguous_assignments[*read.haplotypes->front()].push_back(ambiguous_read_idx);
                    }
                }
                result_.haplotypes[sample].ambiguous_wrt_reference =  result_.haplotypes[sample].ambiguous_wrt_haplotype;
                for (auto& s : possible_ambiguous_assignments) {
                    std::vector<AlignedRead> realigned {};
                    realigned.reserve(s.second.size());
                    for (auto idx : s.second) realigned.push_back(std::move(ambiguous_reads[idx].read));
                    auto& likelihoods = result_.haplotypes[sample].ambiguous_max_likelihoods;
                    safe_realign(realigned, s.first, likelihood_model_, likelihoods, workers);
                    for (std::size_t j {0}; j < s.second.size(); ++j) {
                        result_.haplotypes[sample].ambiguous_wrt_haplotype[s.second[j]].read = realigned[j];
                    }
                    rebase(realigned, s.first, workers);
                    for (std::size_t j {0}; j < s.second.size(); ++j) {
                        result_.haplotypes[sample].ambiguous_wrt_reference[s.second[j]].read = std::move(realigned[j]);
                    }
                }
            }
        }
        for (std::size_t call_idx {0}; call_idx < calls.size(); ++call_idx) {
            std::vector<Allele> alleles {};
            alleles.reserve(calls[call_idx].num_alt() + 1);
            for (auto&& allele : get_resolved_alleles(calls, call_idx, sample)) {
                if (allele) alleles.push_back(std::move(*allele));
            }
            auto allele_support = compute_allele_support(alleles, result_.haplotypes.at(sample), sample, reference);
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

} // namespace csr
} // namespace octopus
