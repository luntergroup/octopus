// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef octopus_vcf_hpp
#define octopus_vcf_hpp

#include <string>

#include <io/variant/vcf_header.hpp>

namespace octopus { namespace vcf {

VcfHeader::Builder make_header_template();

} // namespace vcf
} // namespace octopus

#endif
