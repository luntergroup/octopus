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
                       const Facet::SupportMaps::HaplotypeSupportMaps& support, 
                       const SampleName& sample)
{
    return compute_allele_support(alleles, support.assigned_wrt_reference, support.ambiguous_wrt_reference);
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
    result_.haplotypes.reserve(num_samples);
    for (const auto& p : genotypes) {
        const auto& sample = p.first;
        const auto& sample_genotypes = p.second;
        result_.haplotypes[sample].assigned_wrt_reference.reserve(sample_genotypes.size());
        for (const auto& genotype : sample_genotypes) {
            auto local_reads = copy_overlapped_to_vector(reads.at(sample), genotype);
            for (const auto& haplotype : genotype) {
                // So every called haplotype appears in support map, even if no read support
                result_.haplotypes[sample].assigned_wrt_reference[haplotype] = {};
            }
            if (!local_reads.empty()) {
                // Try to assign each read to a haplotype
                HaplotypeSupportMap genotype_support {};
                if (is_heterozygous(genotype)) {
                    genotype_support = compute_haplotype_support(genotype, local_reads, result_.haplotypes[sample].ambiguous_wrt_haplotype, likelihood_model_);
                } else {
                    if (is_reference(genotype[0])) {
                        genotype_support[genotype[0]] = std::move(local_reads);
                    } else {
                        auto augmented_genotype = genotype;
                        Haplotype ref {mapped_region(genotype), reference};
                        result_.haplotypes[sample].assigned_wrt_reference[ref] = {};
                        augmented_genotype.emplace(std::move(ref));
                        genotype_support = compute_haplotype_support(augmented_genotype, local_reads, result_.haplotypes[sample].ambiguous_wrt_haplotype, likelihood_model_);
                    }
                }
                // Realign assigned reads
                for (auto& s : genotype_support) {
                    const Haplotype& haplotype {s.first};
                    auto& assigned_reads = s.second;
                    safe_realign(assigned_reads, haplotype, likelihood_model_);
                    std::sort(std::begin(assigned_reads), std::end(assigned_reads));
                    result_.haplotypes[sample].assigned_wrt_haplotype[haplotype] = assigned_reads;
                    rebase(assigned_reads, haplotype);
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
                    for (auto idx : s.second) realigned.push_back(ambiguous_reads[idx].read);
                    safe_realign(realigned, s.first, likelihood_model_);
                    for (std::size_t j {0}; j < s.second.size(); ++j) {
                        result_.haplotypes[sample].ambiguous_wrt_haplotype[s.second[j]] = realigned[j];
                    }
                    rebase(realigned, s.first);
                    for (std::size_t j {0}; j < s.second.size(); ++j) {
                        result_.haplotypes[sample].ambiguous_wrt_reference[s.second[j]] = std::move(realigned[j]);
                    }
                }
            }
        }
        for (const auto& call : calls) {
            auto alleles = get_called_alleles(call, sample).first;
            auto allele_support = compute_allele_support(alleles, result_.haplotypes.at(sample), sample);
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
