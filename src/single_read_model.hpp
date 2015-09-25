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
#include <cstddef> // size_t

#include "common.hpp"
#include "haplotype.hpp"
#include "aligned_read.hpp"

namespace Octopus
{
    class SingleReadModel
    {
    public:
        using RealType = Octopus::ProbabilityType;
        
        SingleReadModel()  = delete;
        explicit SingleReadModel(size_t max_cache_size = 10000);
        explicit SingleReadModel(size_t max_num_reads, size_t max_num_haplotypes);
        ~SingleReadModel() = default;
        
        SingleReadModel(const SingleReadModel&)            = default;
        SingleReadModel& operator=(const SingleReadModel&) = default;
        SingleReadModel(SingleReadModel&&)                 = default;
        SingleReadModel& operator=(SingleReadModel&&)      = default;
        
        RealType log_probability(const AlignedRead& read, const Haplotype& haplotype); // ln p(read | haplotype)
        void clear_cache();
        
    private:
        using Cache = std::unordered_map<AlignedRead, std::unordered_map<Haplotype, RealType>>;
        
        Cache cache_;
        size_t max_cache_size_;
        size_t cache_size_;
        
        bool is_cached(const AlignedRead& read, const Haplotype& haplotype) const noexcept;
        void cache(const AlignedRead& read, const Haplotype& haplotype, RealType log_probability);
        RealType get_cached(const AlignedRead& read, const Haplotype& haplotype) const;
    };
    
} // namespace Octopus

#endif /* single_read_model_hpp */
