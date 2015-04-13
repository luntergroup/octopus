//
//  variant_file_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant_file_factory.h"

#include <stdexcept>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "htslib_bcf_facade.h"

namespace fs = boost::filesystem;

std::unique_ptr<IVariantFileImpl> VariantFileFactory::make_impl(const std::string& variant_file_path)
{
    fs::path the_path {variant_file_path};
    
    if (!fs::exists(the_path)) {
        throw std::runtime_error {"Cannot find " + the_path.string()};
    }
    
    return std::make_unique<HtslibBcfFacade>(variant_file_path);
}