//
//  randomiser.hpp
//  Octopus
//
//  Created by Daniel Cooke on 03/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef randomiser_hpp
#define randomiser_hpp

#include <vector>
#include <functional>

#include "candidate_variant_generator.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"

class ReferenceGenome;
class GenomicRegion;

namespace octopus { namespace core { namespace generators
{
    class Randomiser : public CandidateVariantGenerator
    {
    public:
        Randomiser() = delete;
        
        Randomiser(const ReferenceGenome& reference);
        
        Randomiser(const Randomiser&)            = default;
        Randomiser& operator=(const Randomiser&) = default;
        Randomiser(Randomiser&&)                 = default;
        Randomiser& operator=(Randomiser&&)      = default;
        
        ~Randomiser() override = default;
        
        void add_reads(std::vector<AlignedRead>::const_iterator first,
                       std::vector<AlignedRead>::const_iterator last) override;
        void add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                       MappableFlatMultiSet<AlignedRead>::const_iterator last) override;
        
        std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
        
    private:
        std::reference_wrapper<const ReferenceGenome> reference_;
        AlignedRead::RegionType::Size max_read_size_ = 100;
    };
} // namespace generators
} // namespace core
} // namespace octopus

#endif /* randomiser_hpp */
