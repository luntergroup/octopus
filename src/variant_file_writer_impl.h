//
//  variant_file_writer_impl.h
//  Octopus
//
//  Created by Daniel Cooke on 15/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_file_writer_impl__
#define __Octopus__variant_file_writer_impl__

#include <vector>

class Variant;
class GenomicRegion;

class IVariantFileWriterImpl
{
public:
    virtual void write_variants(const std::vector<Variant>& some_variants) = 0;
    virtual ~IVariantFileWriterImpl() = default;
};

#endif /* defined(__Octopus__variant_file_writer_impl__) */
