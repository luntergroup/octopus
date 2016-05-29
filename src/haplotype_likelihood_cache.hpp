//
//  haplotype_likelihood_cache.hpp
//  Octopus
//
//  Created by Daniel Cooke on 25/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef haplotype_likelihood_cache_hpp
#define haplotype_likelihood_cache_hpp

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>

#include "common.hpp"
#include "haplotype.hpp"
#include "aligned_read.hpp"
#include "kmer_mapping.hpp"
#include "haplotype_liklihood_model.hpp"

namespace Octopus
{
    class HaplotypeLikelihoodCache
    {
    public:
        using FlankState = HaplotypeLikelihoodModel::FlankState;
        
        using Likelihoods          = std::vector<double>;
        using LikelihoodsReference = std::reference_wrapper<const Likelihoods>;
        using HaplotypeReference   = std::reference_wrapper<const Haplotype>;
        using SampleLikelihoodMap  = std::unordered_map<HaplotypeReference, LikelihoodsReference>;
        
        HaplotypeLikelihoodCache() = default;
        
        HaplotypeLikelihoodCache(unsigned max_haplotypes, const std::vector<SampleIdType>& samples);
        
        HaplotypeLikelihoodCache(HaplotypeLikelihoodModel likelihood_model,
                                 unsigned max_haplotypes, const std::vector<SampleIdType>& samples);
        
        ~HaplotypeLikelihoodCache() = default;
        
        HaplotypeLikelihoodCache(const HaplotypeLikelihoodCache&)            = default;
        HaplotypeLikelihoodCache& operator=(const HaplotypeLikelihoodCache&) = default;
        HaplotypeLikelihoodCache(HaplotypeLikelihoodCache&&)                 = default;
        HaplotypeLikelihoodCache& operator=(HaplotypeLikelihoodCache&&)      = default;
        
        void populate(const ReadMap& reads, const std::vector<Haplotype>& haplotypes,
                      boost::optional<FlankState> flank_state = boost::none);
        
        const Likelihoods& log_likelihoods(const SampleIdType& sample,
                                           const Haplotype& haplotype) const;
        
        SampleLikelihoodMap extract_sample(const SampleIdType& sample) const;
        
        void reserve(std::size_t num_samples, std::size_t num_haplotypes);
        
        bool contains(const Haplotype& haplotype) const noexcept;
        
        template <typename S, typename Container>
        void insert(S&& sample, const Haplotype& haplotype, Container&& likelihoods);
        
        template <typename Container> void erase(const Container& haplotypes);
        
        void clear();
        
    private:
        static constexpr unsigned char MAPPER_KMER_SIZE {5};
        static constexpr std::size_t MAX_MAPPING_POSITIONS {10};
        
        HaplotypeLikelihoodModel likelihood_model_;
        
        struct ReadPacket
        {
            using Iterator = ReadMap::mapped_type::const_iterator;
            ReadPacket(Iterator first, Iterator last);
            Iterator first, last;
            std::size_t num_reads;
        };
        
        std::unordered_map<Haplotype, std::vector<std::vector<double>>> cache_;
        std::unordered_map<SampleIdType, std::size_t> sample_indices_;
        
        // Just to optimise population
        std::vector<ReadPacket> read_iterators_;
        std::vector<std::size_t> mapping_positions_;
        
        void set_read_iterators_and_sample_indices(const ReadMap& reads);
    };
    
    template <typename S, typename Container>
    void HaplotypeLikelihoodCache::insert(S&& sample, const Haplotype& haplotype,
                                          Container&& likelihoods)
    {
        sample_indices_.emplace(std::forward<S>(sample), sample_indices_.size());
        cache_[haplotype].emplace_back(std::forward<Container>(likelihoods));
    }
    
    template <typename Container>
    void HaplotypeLikelihoodCache::erase(const Container& haplotypes)
    {
        for (const auto& haplotype : haplotypes) {
            cache_.erase(haplotype);
        }
    }
    
    // non-member methods
    
    HaplotypeLikelihoodCache merge_samples(const std::vector<SampleIdType>& samples,
                                           const SampleIdType& new_sample,
                                           const std::vector<Haplotype>& haplotypes,
                                           const HaplotypeLikelihoodCache& haplotype_likelihoods);
    
    namespace debug
    {
        std::vector<std::reference_wrapper<const Haplotype>>
        rank_haplotypes(const std::vector<Haplotype>& haplotypes, const SampleIdType& sample,
                        const HaplotypeLikelihoodCache& haplotype_likelihoods);
        
        template <typename S>
        void print_read_haplotype_liklihoods(S&& stream,
                                             const std::vector<Haplotype>& haplotypes,
                                             const ReadMap& reads,
                                             const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                             const std::size_t n = 5)
        {
            if (n == -1) {
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
            
            using ReadReference = std::reference_wrapper<const AlignedRead>;
            
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
                    ::debug::print_variant_alleles(stream, haplotype);
                    stream << '\n';
                    
                    std::vector<std::pair<ReadReference, double>> likelihoods {};
                    likelihoods.reserve(sample_reads.second.size());
                    
                    std::transform(std::cbegin(sample_reads.second), std::cend(sample_reads.second),
                                   std::cbegin(haplotype_likelihoods.log_likelihoods(sample, haplotype)),
                                   std::back_inserter(likelihoods),
                                   [] (const AlignedRead& read, const double likelihood) {
                                       return std::make_pair(std::cref(read), likelihood);
                                   });
                    
                    const auto mth = std::next(std::begin(likelihoods), m);
                    
                    std::partial_sort(std::begin(likelihoods), mth, std::end(likelihoods),
                                      [] (const auto& lhs, const auto& rhs) {
                                          return lhs.second > rhs.second;
                                      });
                    
                    std::for_each(std::begin(likelihoods), mth,
                                  [&] (const auto& p) {
                                      if (is_single_sample) {
                                          stream << "\t";
                                      } else {
                                          stream << "\t\t";
                                      }
                                      stream << p.first.get().mapped_region()
                                      << " " << p.first.get().cigar_string() << ": ";
                                      stream << p.second << '\n';
                                  });
                }
            }
        }
        
        void print_read_haplotype_liklihoods(const std::vector<Haplotype>& haplotypes,
                                             const ReadMap& reads,
                                             const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                             std::size_t n = 5);
    } // namespace debug
} // namespace Octopus

#endif /* haplotype_likelihood_cache_hpp */
