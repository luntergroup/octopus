// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_filter_utils_hpp
#define variant_filter_utils_hpp

#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"

namespace octopus { namespace csr {

void copy_somatics(const VcfReader& source, VcfWriter& dest);
void copy_somatics(VcfReader::Path source, VcfReader::Path dest);

void copy_denovos(const VcfReader& source, VcfWriter& dest);
void copy_denovos(VcfReader::Path source, VcfReader::Path dest);

} // namespace csr
} // namespace octopus

#endif
