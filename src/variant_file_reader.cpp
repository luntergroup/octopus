//
//  variant_file_reader.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_file_reader.h"

#include "variant.h"
#include "genomic_region.h"

VariantFileReader::VariantFileReader(std::unique_ptr<IVariantFileReaderImpl> the_impl)
:
the_impl_ {std::move(the_impl)}
{}

std::vector<Variant> VariantFileReader::fetch_variants(const GenomicRegion& a_region)
{
    return the_impl_->fetch_variants(a_region);
}
