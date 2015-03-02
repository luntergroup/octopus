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
#include <memory> // std::unique_ptr

#include "variant_file_impl.h"

class Variant;
class GenomicRegion;

class VariantFile
{
public:
    VariantFile() = delete;
    VariantFile(std::unique_ptr<IVariantFileImpl> the_impl);
    ~VariantFile() = default;
    
    VariantFile(const VariantFile&)            = default;
    VariantFile& operator=(const VariantFile&) = default;
    VariantFile(VariantFile&&)                 = default;
    VariantFile& operator=(VariantFile&&)      = default;
    
    std::set<Variant> fetch_variants(const GenomicRegion& a_region);
    void write_variants(const std::set<Variant>& some_variants);
    
private:
    std::unique_ptr<IVariantFileImpl> the_impl_;
};

#endif
