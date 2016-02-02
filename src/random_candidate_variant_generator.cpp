//
//  random_candidate_variant_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "random_candidate_variant_generator.hpp"

#include <random>
#include <chrono>

#include "reference_genome.hpp"
#include "genomic_region.hpp"
#include "variant.hpp"
#include "allele.hpp"
#include "mappable_algorithms.hpp"
#include "sequence_utils.hpp"

namespace Octopus {
    
    RandomCandidateVariantGenerator::RandomCandidateVariantGenerator(const ReferenceGenome& reference)
    :
    reference_ {reference}
    {}
    
    void RandomCandidateVariantGenerator::add_reads(std::vector<AlignedRead>::const_iterator first,
                                                    std::vector<AlignedRead>::const_iterator last)
    {
        max_read_size_ = size(*largest_mappable(first, last));
    }
    
    void RandomCandidateVariantGenerator::add_reads(MappableSet<AlignedRead>::const_iterator first,
                                                    MappableSet<AlignedRead>::const_iterator last)
    {
        max_read_size_ = size(*largest_mappable(first, last));
    }
    
    std::vector<Variant> RandomCandidateVariantGenerator::get_candidates(const GenomicRegion& region)
    {
        auto num_positions = size(region);
        
        std::vector<Variant> result {};
        
        if (num_positions == 0) return result;
        
        static const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        static std::default_random_engine generator {static_cast<unsigned>(seed)};
        
        std::uniform_int_distribution<SizeType> uniform {0, std::min(num_positions, max_read_size_)};
        
        auto positions = decompose(region);
        
        for (auto p = uniform(generator); p < num_positions; p += max_read_size_) {
            auto position = positions[p];
            
            auto reference_allele = get_reference_allele(position, reference_);
            
            Allele mutation {position, reverse_complement_copy(reference_allele.get_sequence())};
            
            result.emplace_back(reference_allele, mutation);
        }
        
        return result;
    }
    
} // namespace Octopus
