// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef haplotype_likelihood_array_hpp
#define haplotype_likelihood_array_hpp

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <iterator>
#include <functional>
#include <iostream>
#include <iomanip>
#include <limits>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/aligned_read.hpp"
#include "basics/aligned_template.hpp"
#include "containers/mappable_block.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "utils/kmer_mapper.hpp"
#include "utils/thread_pool.hpp"
#include "haplotype_likelihood_model.hpp"

namespace octopus {

/*
    HaplotypeLikelihoodArray is essentially a matrix of haplotype likelihoods, i.e.
    p(read | haplotype) for a given set of AlignedReads and Haplotypes.
 
    The matrix can be efficiently populated as the read mapping and alignment are
    done internally which allows minimal memory allocation.
 */
class HaplotypeLikelihoodArray
{
public:
    using FlankState = HaplotypeLikelihoodModel::FlankState;
    
    using LogProbability       = HaplotypeLikelihoodModel::LogProbability;
    using LikelihoodVector     = std::vector<LogProbability>;
    using LikelihoodVectorRef  = std::reference_wrapper<const LikelihoodVector>;
    using HaplotypeRef         = std::reference_wrapper<const Haplotype>;
    using SampleLikelihoodMap  = std::unordered_map<HaplotypeRef, LikelihoodVectorRef>;

    using OptionalThreadPool = boost::optional<ThreadPool&>;
    
    HaplotypeLikelihoodArray() = default;
    
    HaplotypeLikelihoodArray(unsigned num_haplotypes_hint, const std::vector<SampleName>& samples);
    
    HaplotypeLikelihoodArray(HaplotypeLikelihoodModel likelihood_model,
                             unsigned num_haplotypes_hint,
                             const std::vector<SampleName>& samples);
    
    HaplotypeLikelihoodArray(const HaplotypeLikelihoodArray&)            = default;
    HaplotypeLikelihoodArray& operator=(const HaplotypeLikelihoodArray&) = default;
    HaplotypeLikelihoodArray(HaplotypeLikelihoodArray&&)                 = default;
    HaplotypeLikelihoodArray& operator=(HaplotypeLikelihoodArray&&)      = default;
    
    ~HaplotypeLikelihoodArray() = default;
    
    void populate(const ReadMap& reads,
                  const MappableBlock<Haplotype>& haplotypes,
                  boost::optional<FlankState> flank_state = boost::none,
                  OptionalThreadPool workers = boost::none);
    void populate(const TemplateMap& reads,
                  const MappableBlock<Haplotype>& haplotypes,
                  boost::optional<FlankState> flank_state = boost::none,
                  OptionalThreadPool workers = boost::none);
    
    std::size_t num_likelihoods(const SampleName& sample) const;
    std::size_t num_likelihoods() const; // if prmed
    
    const LikelihoodVector& operator()(const SampleName& sample, const Haplotype& haplotype) const;
    const LikelihoodVector& operator()(const SampleName& sample, const IndexedHaplotype<>& haplotype) const;
    const LikelihoodVector& operator[](const Haplotype& haplotype) const; // when primed with a sample
    const LikelihoodVector& operator[](const IndexedHaplotype<>& haplotype) const noexcept; // when primed with a sample
    
    std::vector<SampleName> samples() const;
    MappableBlock<Haplotype> haplotypes() const;
    
    SampleLikelihoodMap extract_sample(const SampleName& sample) const;
    
    bool contains(const Haplotype& haplotype) const noexcept;
    
    void reset(MappableBlock<Haplotype> haplotypes);
    
    bool is_empty() const noexcept;
    
    void clear() noexcept;
    
    bool is_primed() const noexcept;
    void prime(const SampleName& sample) const;
    void unprime() const noexcept;
    
    HaplotypeLikelihoodArray merge_samples(const std::vector<SampleName>& samples, boost::optional<SampleName> new_sample = boost::none) const;
    HaplotypeLikelihoodArray merge_samples(boost::optional<SampleName> new_sample = boost::none) const;
    
private:
    static constexpr unsigned char mapperKmerSize {6};
    static constexpr std::size_t maxMappingPositions {10};
    
    HaplotypeLikelihoodModel likelihood_model_;
    
    struct ReadPacket
    {
        using Iterator = ReadMap::mapped_type::const_iterator;
        ReadPacket(Iterator first, Iterator last);
        Iterator first, last;
        std::size_t num_reads;
    };
    struct TemplatePacket
    {
        using Iterator = TemplateMap::mapped_type::const_iterator;
        TemplatePacket(Iterator first, Iterator last);
        Iterator first, last;
        std::size_t num_templates;
    };
    
    std::vector<std::vector<LikelihoodVector>> likelihoods_;
    std::unordered_map<Haplotype, std::size_t, HaplotypeHash> haplotype_indices_;
    std::unordered_map<SampleName, std::size_t> sample_indices_;
    std::vector<SampleName> samples_;
    MappableBlock<Haplotype> haplotypes_;
    
    mutable boost::optional<std::size_t> primed_sample_;
    
    // Just to optimise population
    std::vector<ReadPacket> read_iterators_;
    std::vector<TemplatePacket> template_iterators_;
    std::vector<std::size_t> mapping_positions_;
    
    void set_read_iterators_and_sample_indices(const ReadMap& reads);
    void set_template_iterators_and_sample_indices(const TemplateMap& reads);
};

// non-member methods

namespace debug {

std::vector<std::reference_wrapper<const Haplotype>>
rank_haplotypes(const MappableBlock<Haplotype>& haplotypes, const SampleName& sample,
                const HaplotypeLikelihoodArray& haplotype_likelihoods);

template <typename S>
void print_read_likelihood_helper(S& stream, const AlignedRead& read, const double likelihood)
{
    stream << read.name() << " "
           << mapped_region(read) << " "
           << read.cigar() << ": "
           << likelihood << '\n';
}
template <typename S>
void print_read_likelihood_helper(S& stream, const AlignedTemplate& reads, const double likelihood)
{
    stream << "[ ";
    for (const auto& read : reads) {
        stream <<  "<" << read.name() << " " << mapped_region(read) << " " << read.cigar() << "> "; 
    }
    stream << "] : " << likelihood << '\n';
}

template <typename S, typename Read>
void print_read_haplotype_likelihoods(S&& stream,
                                     const MappableBlock<Haplotype>& haplotypes,
                                     const MappableMap<SampleName, Read>& reads,
                                     const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                     const std::size_t n = std::numeric_limits<std::size_t>::max())
{
    if (n == static_cast<std::size_t>(-1)) {
        stream << "Printing all read likelihoods for each haplotype in ";
    } else {
        stream << "Printing top " << n << " read likelihoods for each haplotype in ";
    }
    const bool is_single_sample {reads.size() == 1};
    if (is_single_sample) {
        stream << "sample " << std::cbegin(reads)->first;
    } else {
        stream << "each sample";
    }
    stream << '\n';
    using ReadReference = std::reference_wrapper<const Read>;
    for (const auto& sample_reads : reads) {
        const auto& sample = sample_reads.first;
        if (!is_single_sample) {
            stream << "Sample: " << sample << ":" << '\n';
        }
        const auto ranked_haplotypes = rank_haplotypes(haplotypes, sample, haplotype_likelihoods);
        const auto m = std::min(n, sample_reads.second.size());
        for (const auto& haplotype : ranked_haplotypes) {
            if (!is_single_sample) {
                stream << "\t";
            }
            debug::print_variant_alleles(stream, haplotype);
            stream << '\n';
            std::vector<std::pair<ReadReference, HaplotypeLikelihoodArray::LogProbability >> likelihoods {};
            likelihoods.reserve(sample_reads.second.size());
            std::transform(std::cbegin(sample_reads.second), std::cend(sample_reads.second),
                           std::cbegin(haplotype_likelihoods(sample, haplotype)),
                           std::back_inserter(likelihoods),
                           [] (const Read& read, const auto likelihood) {
                               return std::make_pair(std::cref(read), likelihood);
                           });
            const auto mth = std::next(std::begin(likelihoods), m);
            std::partial_sort(std::begin(likelihoods), mth, std::end(likelihoods),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.second > rhs.second;
                              });
            std::for_each(std::begin(likelihoods), mth, [&stream, is_single_sample] (const auto& p) {
                if (is_single_sample) {
                    stream << "\t";
                } else {
                    stream << "\t\t";
                }
                print_read_likelihood_helper(stream, p.first.get(), p.second);
            });
        }
    }
}

template <typename Read>
void print_read_haplotype_likelihoods(const MappableBlock<Haplotype>& haplotypes,
                                      const MappableMap<SampleName, Read>& reads,
                                      const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                      std::size_t n = std::numeric_limits<std::size_t>::max())
{
    print_read_haplotype_likelihoods(std::cout, haplotypes, reads, haplotype_likelihoods, n);
}

} // namespace debug

} // namespace octopus

#endif
