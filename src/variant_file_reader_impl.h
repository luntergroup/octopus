//
//  variant_file_reader_impl.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_file_reader_impl_h
#define Octopus_variant_file_reader_impl_h

#include <vector>

class Variant;
class GenomicRegion;

class IVariantFileReaderImpl
{
public:
    virtual std::vector<Variant> fetch_variants(const GenomicRegion& a_region) = 0;
    virtual ~IVariantFileReaderImpl() = default;
};

#endif
