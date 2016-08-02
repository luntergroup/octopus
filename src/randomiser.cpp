//
//  randomiser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "randomiser.hpp"

#include <random>
#include <chrono>

#include "reference_genome.hpp"
#include "genomic_region.hpp"
#include "allele.hpp"
#include "mappable_algorithms.hpp"
#include "sequence_utils.hpp"

namespace octopus { namespace coretools {

Randomiser::Randomiser(const ReferenceGenome& reference)
:
reference_ {reference}
{}

std::unique_ptr<VariantGenerator> Randomiser::do_clone() const
{
    return std::make_unique<Randomiser>(*this);
}

void Randomiser::do_add_reads(VectorIterator first, VectorIterator last)
{
    max_read_size_ = region_size(*largest_mappable(first, last));
}

void Randomiser::do_add_reads(FlatSetIterator first, FlatSetIterator last)
{
    max_read_size_ = region_size(*largest_mappable(first, last));
}

std::vector<Variant> Randomiser::do_generate_variants(const GenomicRegion& region)
{
    auto num_positions = region_size(region);
    
    std::vector<Variant> result {};
    
    if (num_positions == 0) return result;
    
    static const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    
    static std::default_random_engine generator {static_cast<unsigned>(seed)};
    
    using T = Variant::RegionType::Size;
    
    std::uniform_int_distribution<T> uniform {0, std::min(num_positions, max_read_size_)};
    
    auto positions = decompose(region);
    
    for (auto p = uniform(generator); p < num_positions; p += max_read_size_) {
        auto position = positions[p];
        
        auto reference_allele = make_reference_allele(position, reference_);
        
        Allele mutation {position, utils::reverse_complement_copy(reference_allele.sequence())};
        
        result.emplace_back(reference_allele, mutation);
    }
    
    return result;
}

std::string Randomiser::name() const
{
    return "Random";
}

} // namespace coretools
} // namespace octopus
