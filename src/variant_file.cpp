//
//  variant_file.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_file.h"

#include "variant.h"
#include "genomic_region.h"

VariantFile::VariantFile(std::unique_ptr<IVariantFileImpl> the_impl)
: the_impl_ {std::move(the_impl)}
{}

std::vector<Variant> VariantFile::fetch_variants(const GenomicRegion& a_region)
{
    return the_impl_->fetch_variants(a_region);
}

void VariantFile::write_variants(const std::vector<Variant>& some_variants)
{
    return the_impl_->write_variants(some_variants);
}
