//
//  variant_file_implementer.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_file_impl_h
#define Octopus_variant_file_impl_h

#include <vector>

class Variant;
class GenomicRegion;

class IVariantFileImpl
{
public:
    virtual std::vector<Variant> fetch_variants(const GenomicRegion& a_region) = 0;
    virtual void write_variants(const std::vector<Variant>& some_variants) = 0;
    virtual ~IVariantFileImpl() = default;
};

#endif
