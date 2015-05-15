//
//  variant_file_reader.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_file_reader_h
#define Octopus_variant_file_reader_h

#include <vector>
#include <memory> // std::unique_ptr

#include "variant_file_reader_impl.h"

class Variant;
class GenomicRegion;

class VariantFileReader
{
public:
    VariantFileReader() = delete;
    VariantFileReader(std::unique_ptr<IVariantFileReaderImpl> the_impl);
    ~VariantFileReader() = default;
    
    VariantFileReader(const VariantFileReader&)            = default;
    VariantFileReader& operator=(const VariantFileReader&) = default;
    VariantFileReader(VariantFileReader&&)                 = default;
    VariantFileReader& operator=(VariantFileReader&&)      = default;
    
    std::vector<Variant> fetch_variants(const GenomicRegion& a_region);
    
private:
    std::unique_ptr<IVariantFileReaderImpl> the_impl_;
};

#endif
