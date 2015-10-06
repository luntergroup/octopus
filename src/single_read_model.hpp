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
        SingleReadModel()  = delete;
        explicit SingleReadModel(size_t max_num_reads, size_t max_num_haplotypes);
        ~SingleReadModel() = default;
        
        SingleReadModel(const SingleReadModel&)            = default;
        SingleReadModel& operator=(const SingleReadModel&) = default;
        SingleReadModel(SingleReadModel&&)                 = default;
        SingleReadModel& operator=(SingleReadModel&&)      = default;
        
        double log_probability(const AlignedRead& read, const Haplotype& haplotype); // ln p(read | haplotype)
        
        void clear_cache();
        
    private:
        using Cache = std::unordered_map<AlignedRead, std::unordered_map<Haplotype, double>>;
        
        Cache cache_;
        size_t max_num_reads_;
        size_t max_num_haplotypes_;
        
        bool is_cached(const AlignedRead& read, const Haplotype& haplotype) const noexcept;
        void cache(const AlignedRead& read, const Haplotype& haplotype, double value);
        double get_cached(const AlignedRead& read, const Haplotype& haplotype) const;
    };
    
} // namespace Octopus

#endif /* single_read_model_hpp */
