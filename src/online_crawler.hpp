//
//  online_crawler.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/06/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__online_crawler__
#define __Octopus__online_crawler__

#include <vector>
#include <functional>

#include "candidate_variant_generator.hpp"
#include "variant.hpp"

class GenomicRegion;
class ReferenceGenome;
class AlignedRead;

namespace octopus { namespace core { namespace generators
{
class OnlineCrawler : public CandidateVariantGenerator
{
public:
    OnlineCrawler() = delete;
    
    OnlineCrawler(const ReferenceGenome& reference, Variant::RegionType::Size max_variant_size = 100);
    
    OnlineCrawler(const OnlineCrawler&)            = default;
    OnlineCrawler& operator=(const OnlineCrawler&) = default;
    OnlineCrawler(OnlineCrawler&&)                 = default;
    OnlineCrawler& operator=(OnlineCrawler&&)      = default;
    
    ~OnlineCrawler() override = default;
    
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;

private:
    std::reference_wrapper<const ReferenceGenome> reference_;
    Variant::RegionType::Size max_variant_size_;
};
} // namespace generators
} // namespace core
} // namespace octopus

#endif /* defined(__Octopus__online_crawler__) */
