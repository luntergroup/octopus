//
//  variant_file.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_file_h
#define Octopus_variant_file_h

#include <set>

#include "variant.h"
#include "genomic_region.h"

class VariantFile
{
public:
    std::set<Variant> fetch_variants(const GenomicRegion& a_region);
private:
};

inline
std::set<Variant> VariantFile::fetch_variants(const GenomicRegion& a_region)
{
    return std::set<Variant> {};
}

#endif
