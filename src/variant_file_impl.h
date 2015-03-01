//
//  variant_file_implementer.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_file_impl_h
#define Octopus_variant_file_impl_h

#include <set>

class Variant;
class GenomicRegion;

class VariantFileImpl
{
public:
    virtual std::set<Variant> fetch_variants(const GenomicRegion& a_region);
};

#endif
