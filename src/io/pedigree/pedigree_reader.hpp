// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef pedigree_reader_hpp
#define pedigree_reader_hpp

#include <boost/filesystem/path.hpp>

#include "basics/pedigree.hpp"

namespace octopus {namespace io {

Pedigree read_pedigree(const boost::filesystem::path& ped_file);

} // namespace io
} // namespace octopus

#endif
