//
//  variant_caller.h
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_caller_h
#define Octopus_variant_caller_h

#include <vector>

class AlignedRead;
class VcfRecord;

class VariantCaller
{
public:
    std::vector<VcfRecord> call_variants(const std::vector<AlignedRead>& reads);
    
    virtual ~VariantCaller() = default;
    
private:
    
};

#endif
