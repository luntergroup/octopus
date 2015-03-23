//
//  haplotype.h
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype__
#define __Octopus__haplotype__

#include <queue>

#include "variant.h"

class ReferenceGenome;

using Haplotype = std::deque<Variant>;

//class Haplotype
//{
//public:
//    Haplotype() = delete;
//    explicit Haplotype(ReferenceGenome& the_reference);
//    ~Haplotype() = default;
//    
//    Haplotype(const Haplotype&)            = default;
//    Haplotype& operator=(const Haplotype&) = default;
//    Haplotype(Haplotype&&)                 = default;
//    Haplotype& operator=(Haplotype&&)      = default;
//    
//    
//private:
//    ReferenceGenome& the_reference_;
//};

#endif /* defined(__Octopus__haplotype__) */
