//
//  single_read_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef single_read_model_hpp
#define single_read_model_hpp

#include <unordered_map>
#include <functional>
#include <cstddef> // size_t

#include "common.hpp"
#include "haplotype.hpp"
#include "aligned_read.hpp"

namespace std
{
    template <typename T> struct hash<reference_wrapper<const T>>
    {
        size_t operator()(reference_wrapper<const T> r) const
        {
            return hash<T>()(r);
        }
    };
} // namespace std

namespace Octopus
{
    class SingleReadModel
    {
    public:
        SingleReadModel()  = default;
        explicit SingleReadModel(size_t max_num_reads, size_t max_num_haplotypes);
        ~SingleReadModel() = default;
        
        SingleReadModel(const SingleReadModel&)            = default;
        SingleReadModel& operator=(const SingleReadModel&) = default;
        SingleReadModel(SingleReadModel&&)                 = default;
        SingleReadModel& operator=(SingleReadModel&&)      = default;
        
        double log_probability(const AlignedRead& read, const Haplotype& haplotype); // ln p(read | haplotype)
        
        void clear_cache();
        
    private:
        using ReadReference      = std::reference_wrapper<const AlignedRead>;
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        using Cache = std::unordered_map<ReadReference, std::unordered_map<HaplotypeReference, double>>;
        
        Cache cache_;
        size_t max_num_reads_;
        size_t max_num_haplotypes_;
        
        bool is_cached(const AlignedRead& read, const Haplotype& haplotype) const noexcept;
        void cache(const AlignedRead& read, const Haplotype& haplotype, double value);
        double get_cached(const AlignedRead& read, const Haplotype& haplotype) const;
    };
    
} // namespace Octopus

#endif /* single_read_model_hpp */
