//
//  variant_file_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_file_factory__
#define __Octopus__variant_file_factory__

#include <string>
#include <memory> // std::unique_ptr

#include "variant_file_reader_impl.h"

class VariantFileFactory
{
public:
    VariantFileFactory() = default;
    
    std::unique_ptr<IVariantFileReaderImpl> make_reader(const std::string& variant_file_path);
};

#endif /* defined(__Octopus__variant_file_factory__) */
