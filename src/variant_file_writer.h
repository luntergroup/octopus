//
//  variant_file_writer.h
//  Octopus
//
//  Created by Daniel Cooke on 15/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_file_writer__
#define __Octopus__variant_file_writer__

#include <vector>
#include <memory> // std::unique_ptr

#include "variant_file_writer_impl.h"

class Variant;
class GenomicRegion;

class VariantFileWriter
{
public:
    VariantFileWriter() = delete;
    VariantFileWriter(std::unique_ptr<IVariantFileWriterImpl> the_impl);
    ~VariantFileWriter() = default;
    
    VariantFileWriter(const VariantFileWriter&)            = default;
    VariantFileWriter& operator=(const VariantFileWriter&) = default;
    VariantFileWriter(VariantFileWriter&&)                 = default;
    VariantFileWriter& operator=(VariantFileWriter&&)      = default;
    
    void write_variants(const std::vector<Variant>& some_variants);
    
private:
    std::unique_ptr<IVariantFileWriterImpl> the_impl_;
};

#endif /* defined(__Octopus__variant_file_writer__) */
