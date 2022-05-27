// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "haplotype_likelihood_array.hpp"

#include <utility>
#include <cassert>
#include <deque>

#include "utils/erase_if.hpp"

namespace octopus {

// public methods

HaplotypeLikelihoodArray::HaplotypeLikelihoodArray(const unsigned num_haplotypes_hint,
                                                   const std::vector<SampleName>& samples)
: likelihoods_ {}
, haplotype_indices_ {num_haplotypes_hint}
, sample_indices_ {samples.size()}
, samples_ {samples}
{
    mapping_positions_.resize(maxMappingPositions);
}

HaplotypeLikelihoodArray::HaplotypeLikelihoodArray(HaplotypeLikelihoodModel likelihood_model,
                                                   unsigned num_haplotypes_hint,
                                                   const std::vector<SampleName>& samples)
: likelihood_model_ {std::move(likelihood_model)}
, likelihoods_ {}
, haplotype_indices_ {num_haplotypes_hint}
, sample_indices_ {samples.size()}
, samples_ {samples}
{
    mapping_positions_.resize(maxMappingPositions);
}

HaplotypeLikelihoodArray::ReadPacket::ReadPacket(Iterator first, Iterator last)
: first {first}
, last {last}
, num_reads {static_cast<std::size_t>(std::distance(first, last))}
{}

HaplotypeLikelihoodArray::TemplatePacket::TemplatePacket(Iterator first, Iterator last)
: first {first}
, last {last}
, num_templates {static_cast<std::size_t>(std::distance(first, last))}
{}

void HaplotypeLikelihoodArray::populate(const ReadMap& reads,
                                        const MappableBlock<Haplotype>& haplotypes,
                                        boost::optional<FlankState> flank_state,
                                        OptionalThreadPool workers)
{
    // This code is not very pretty because it is a bottleneck for the entire application.
    // We want to try a minimise memory allocations for the mapping.
    haplotype_indices_.clear();
    if (haplotype_indices_.bucket_count() < haplotypes.size()) {
        haplotype_indices_.rehash(haplotypes.size());
    }
    set_read_iterators_and_sample_indices(reads);
    assert(reads.size() == read_iterators_.size());
    const auto num_samples = reads.size();
    // Precompute all read hashes so we don't have to recompute for each haplotype
    std::vector<std::vector<KmerPerfectHashes>> read_hashes {};
    read_hashes.reserve(num_samples);
    for (const auto& t : read_iterators_) {
        std::vector<KmerPerfectHashes> sample_read_hashes {};
        sample_read_hashes.reserve(t.num_reads);
        std::transform(t.first, t.last, std::back_inserter(sample_read_hashes),
                       [] (const AlignedRead& read) { return compute_kmer_hashes<mapperKmerSize>(read.sequence()); });
        read_hashes.emplace_back(std::move(sample_read_hashes));
    }
    auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
    const auto first_mapping_position = std::begin(mapping_positions_);
    likelihoods_.resize(haplotypes.size(), std::vector<LikelihoodVector>(num_samples));
    for (std::size_t haplotype_idx {0}; haplotype_idx < haplotypes.size(); ++haplotype_idx) {
        const auto& haplotype = haplotypes[haplotype_idx];
        populate_kmer_hash_table<mapperKmerSize>(haplotype.sequence(), haplotype_hashes);
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        likelihood_model_.reset(haplotype, flank_state);
        for (std::size_t sample_idx {0}; sample_idx < num_samples; ++sample_idx) {
            auto& likelihoods = likelihoods_[haplotype_idx][sample_idx];
            const auto& t = read_iterators_[sample_idx];
            likelihoods.resize(t.num_reads);
            std::transform(t.first, t.last, std::cbegin(read_hashes[sample_idx]), std::begin(likelihoods),
                           [&] (const AlignedRead& read, const auto& read_hashes) {
                               const auto last_mapping_position = map_query_to_target(read_hashes, haplotype_hashes,
                                                                                      haplotype_mapping_counts,
                                                                                      first_mapping_position,
                                                                                      maxMappingPositions);
                               reset_mapping_counts(haplotype_mapping_counts);
                               return likelihood_model_.evaluate(read, first_mapping_position, last_mapping_position);
                           });
        }
        clear_kmer_hash_table(haplotype_hashes);
        haplotype_indices_.emplace(haplotype, haplotype_idx);
    }
    likelihood_model_.clear();
    read_iterators_.clear();
    haplotypes_ = haplotypes;
}

void HaplotypeLikelihoodArray::populate(const TemplateMap& reads,
                                        const MappableBlock<Haplotype>& haplotypes,
                                        boost::optional<FlankState> flank_state,
                                        OptionalThreadPool workers)
{
    haplotype_indices_.clear();
    if (haplotype_indices_.bucket_count() < haplotypes.size()) {
        haplotype_indices_.rehash(haplotypes.size());
    }
    set_template_iterators_and_sample_indices(reads);
    assert(reads.size() == template_iterators_.size());
    const auto num_samples = reads.size();
    // Precompute all read hashes so we don't have to recompute for each haplotype
    std::vector<std::vector<std::vector<KmerPerfectHashes>>> template_hashes {};
    template_hashes.reserve(num_samples);
    for (const auto& t : template_iterators_) {
        std::vector<std::vector<KmerPerfectHashes>> sample_read_hashes {};
        sample_read_hashes.reserve(t.num_templates);
        std::transform(t.first, t.last, std::back_inserter(sample_read_hashes), [] (const AlignedTemplate& reads) {
            std::vector<KmerPerfectHashes> result {};
            result.reserve(reads.size());
            for (const auto& read : reads) result.push_back(compute_kmer_hashes<mapperKmerSize>(read.sequence()));
            return result;
        });
        template_hashes.emplace_back(std::move(sample_read_hashes));
    }
    auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
    thread_local std::vector<HaplotypeLikelihoodModel::MappingPositionVector> mapping_positions {};
    likelihoods_.resize(haplotypes.size(), std::vector<LikelihoodVector>(num_samples));
    for (std::size_t haplotype_idx {0}; haplotype_idx < haplotypes.size(); ++haplotype_idx) {
        const auto& haplotype = haplotypes[haplotype_idx];
        populate_kmer_hash_table<mapperKmerSize>(haplotype.sequence(), haplotype_hashes);
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        likelihood_model_.reset(haplotype, flank_state);
        for (std::size_t sample_idx {0}; sample_idx < num_samples; ++sample_idx) {
            auto& likelihoods = likelihoods_[haplotype_idx][sample_idx];
            const auto& t = template_iterators_[sample_idx];
            likelihoods.resize(t.num_templates);
            std::transform(t.first, t.last, std::cbegin(template_hashes[sample_idx]), std::begin(likelihoods),
                           [&] (const AlignedTemplate& read_template, const auto& template_hashes) {
                               mapping_positions.resize(read_template.size());
                               assert(read_template.size() == template_hashes.size());
                               for (std::size_t i {0}; i < template_hashes.size(); ++i) {
                                   mapping_positions[i].resize(maxMappingPositions);
                                   mapping_positions[i].erase(map_query_to_target(template_hashes[i], haplotype_hashes,
                                                                                  haplotype_mapping_counts,
                                                                                  std::begin(mapping_positions[i]),
                                                                                  maxMappingPositions),
                                                              std::end(mapping_positions[i]));
                                   reset_mapping_counts(haplotype_mapping_counts);
                               }
                               return likelihood_model_.evaluate(read_template, mapping_positions);
                           });
        }
        clear_kmer_hash_table(haplotype_hashes);
        haplotype_indices_.emplace(haplotype, haplotype_idx);
    }
    likelihood_model_.clear();
    read_iterators_.clear();
    haplotypes_ = haplotypes;
}

std::size_t HaplotypeLikelihoodArray::num_likelihoods(const SampleName& sample) const
{
    return likelihoods_.front()[sample_indices_.at(sample)].size();
}

std::size_t HaplotypeLikelihoodArray::num_likelihoods() const // if primed
{
    assert(is_primed());
    return likelihoods_.front()[*primed_sample_].size();
}

const HaplotypeLikelihoodArray::LikelihoodVector&
HaplotypeLikelihoodArray::operator()(const SampleName& sample, const Haplotype& haplotype) const
{
    return likelihoods_[haplotype_indices_.at(haplotype)][sample_indices_.at(sample)];
}

const HaplotypeLikelihoodArray::LikelihoodVector&
HaplotypeLikelihoodArray::operator()(const SampleName& sample, const IndexedHaplotype<>& haplotype) const
{
    return likelihoods_[index_of(haplotype)][sample_indices_.at(sample)];
}

const HaplotypeLikelihoodArray::LikelihoodVector&
HaplotypeLikelihoodArray::operator[](const Haplotype& haplotype) const
{
    assert(is_primed());
    return likelihoods_[haplotype_indices_.at(haplotype)][*primed_sample_];
}

const HaplotypeLikelihoodArray::LikelihoodVector&
HaplotypeLikelihoodArray::operator[](const IndexedHaplotype<>& haplotype) const noexcept
{
    assert(is_primed());
    return likelihoods_[index_of(haplotype)][*primed_sample_];
}

std::vector<SampleName> HaplotypeLikelihoodArray::samples() const
{
    return samples_;
}

MappableBlock<Haplotype> HaplotypeLikelihoodArray::haplotypes() const
{
    return haplotypes_;
}

HaplotypeLikelihoodArray::SampleLikelihoodMap
HaplotypeLikelihoodArray::extract_sample(const SampleName& sample) const
{
    const auto sample_index = sample_indices_.at(sample);
    SampleLikelihoodMap result {haplotype_indices_.size()};
    for (const auto& p : haplotype_indices_) {
        result.emplace(p.first, likelihoods_[p.second][sample_index]);
    }
    return result;
}

bool HaplotypeLikelihoodArray::contains(const Haplotype& haplotype) const noexcept
{
    return haplotype_indices_.count(haplotype) == 1;
}

bool HaplotypeLikelihoodArray::is_empty() const noexcept
{
    return likelihoods_.empty();
}

void HaplotypeLikelihoodArray::clear() noexcept
{
    likelihoods_.clear();
    haplotype_indices_.clear();
    sample_indices_.clear();
    haplotypes_.clear();
    unprime();
}

bool HaplotypeLikelihoodArray::is_primed() const noexcept
{
    return static_cast<bool>(primed_sample_);
}

void HaplotypeLikelihoodArray::prime(const SampleName& sample) const
{
    primed_sample_ = sample_indices_.at(sample);
}

void HaplotypeLikelihoodArray::unprime() const noexcept
{
    primed_sample_ = boost::none;
}

// private methods

void HaplotypeLikelihoodArray::set_read_iterators_and_sample_indices(const ReadMap& reads)
{
    read_iterators_.clear();
    sample_indices_.clear();
    const auto num_samples = reads.size();
    if (read_iterators_.capacity() < num_samples) {
        read_iterators_.reserve(num_samples);
    }
    if (sample_indices_.bucket_count() < num_samples) {
        sample_indices_.rehash(num_samples);
    }
    std::size_t i {0};
    for (const auto& p : reads) {
        read_iterators_.emplace_back(std::cbegin(p.second), std::cend(p.second));
        sample_indices_.emplace(p.first, i++);
    }
}

void HaplotypeLikelihoodArray::set_template_iterators_and_sample_indices(const TemplateMap& reads)
{
    template_iterators_.clear();
    sample_indices_.clear();
    const auto num_samples = reads.size();
    if (template_iterators_.capacity() < num_samples) {
        template_iterators_.reserve(num_samples);
    }
    if (sample_indices_.bucket_count() < num_samples) {
        sample_indices_.rehash(num_samples);
    }
    std::size_t i {0};
    for (const auto& p : reads) {
        template_iterators_.emplace_back(std::cbegin(p.second), std::cend(p.second));
        sample_indices_.emplace(p.first, i++);
    }
}

void HaplotypeLikelihoodArray::reset(MappableBlock<Haplotype> haplotypes)
{
    assert(haplotypes.size() <= haplotypes_.size());
    if (haplotypes.empty()) {
        clear();
    } else if (haplotypes.size() < haplotypes_.size()) {
        std::deque<std::size_t> indices_to_keep {};
        for (std::size_t haplotype_idx {0}; haplotype_idx < haplotypes.size(); ++haplotype_idx) {
            auto& old_haplotype_idx = haplotype_indices_.at(haplotypes[haplotype_idx]);
            indices_to_keep.push_back(old_haplotype_idx);
            old_haplotype_idx = haplotype_idx;
        }
        assert(!indices_to_keep.empty());
        assert(std::is_sorted(std::cbegin(indices_to_keep), std::cend(indices_to_keep)));
        std::size_t idx {0};
        const auto remover = [&] (const auto&) noexcept {
            if (indices_to_keep.empty()) return true;
            if (idx++ != indices_to_keep.front()) return true;
            indices_to_keep.pop_front();
            return false;
        };
        erase_if(likelihoods_, remover);
        haplotypes_ = std::move(haplotypes);
    }
}

HaplotypeLikelihoodArray
HaplotypeLikelihoodArray::merge_samples(const std::vector<SampleName>& samples, boost::optional<SampleName> new_sample) const
{
    if (!new_sample) {
        new_sample = SampleName {};
        for (const auto& sample : samples) *new_sample += sample;
    }
    HaplotypeLikelihoodArray result {static_cast<unsigned>(haplotypes_.size()), {std::move(*new_sample)}};
    result.haplotypes_ = haplotypes_;
    std::size_t total_num_likelihoods {0};
    for (const auto& sample : samples) {
        total_num_likelihoods += this->num_likelihoods(sample);
    }
    result.likelihoods_.resize(haplotypes_.size(), std::vector<LikelihoodVector>(1, LikelihoodVector(total_num_likelihoods)));
    for (std::size_t haplotype_idx {0}; haplotype_idx < haplotypes_.size(); ++haplotype_idx) {
        auto dst_likelihood_itr = std::begin(result.likelihoods_[haplotype_idx][0]);
        for (const auto& sample : samples) {
            const auto& src_likelihoods = likelihoods_[haplotype_idx][sample_indices_.at(sample)];
            dst_likelihood_itr = std::copy(std::cbegin(src_likelihoods), std::cend(src_likelihoods), dst_likelihood_itr);
        }
    }
    result.haplotype_indices_ = haplotype_indices_;
    result.sample_indices_.emplace(result.samples_.front(), 0);
    result.prime(result.samples_.front());
    return result;
}

HaplotypeLikelihoodArray
HaplotypeLikelihoodArray::merge_samples(boost::optional<SampleName> new_sample) const
{
    if (!new_sample) {
        new_sample = SampleName {};
        for (const auto& sample : samples_) *new_sample += sample;
    }
    HaplotypeLikelihoodArray result {static_cast<unsigned>(haplotypes_.size()), {std::move(*new_sample)}};
    result.haplotypes_ = haplotypes_;
    std::size_t total_num_likelihoods {0};
    for (const auto& s : likelihoods_.front()) {
        total_num_likelihoods += s.size();
    }
    result.likelihoods_.resize(haplotypes_.size(), std::vector<LikelihoodVector>(1, LikelihoodVector(total_num_likelihoods)));
    for (std::size_t haplotype_idx {0}; haplotype_idx < haplotypes_.size(); ++haplotype_idx) {
        auto dst_likelihood_itr = std::begin(result.likelihoods_[haplotype_idx][0]);
        for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
            const auto& src_likelihoods = likelihoods_[haplotype_idx][sample_idx];
            dst_likelihood_itr = std::copy(std::cbegin(src_likelihoods), std::cend(src_likelihoods), dst_likelihood_itr);
        }
    }
    result.haplotype_indices_ = haplotype_indices_;
    result.sample_indices_.emplace(result.samples_.front(), 0);
    result.prime(result.samples_.front());
    return result;
}

// non-member methods

namespace debug {

std::vector<std::reference_wrapper<const Haplotype>>
rank_haplotypes(const MappableBlock<Haplotype>& haplotypes,
                const SampleName& sample,
                const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    std::vector<std::pair<std::reference_wrapper<const Haplotype>, HaplotypeLikelihoodArray::LogProbability>> ranks {};
    ranks.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        const auto& likelihoods = haplotype_likelihoods(sample, haplotype);
        ranks.emplace_back(haplotype, std::accumulate(std::cbegin(likelihoods), std::cend(likelihoods), 0.0));
    }
    std::sort(std::begin(ranks), std::end(ranks),
              [] (const auto& lhs, const auto& rhs) {
                  return lhs.second > rhs.second;
              });
    std::vector<std::reference_wrapper<const Haplotype>> result {};
    result.reserve(haplotypes.size());
    std::transform(std::cbegin(ranks), std::cend(ranks), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    return result;
}

} // namespace debug

} // namespace octopus
