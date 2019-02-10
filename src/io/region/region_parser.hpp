// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef region_parser_hpp
#define region_parser_hpp

#include <string>
#include <deque>

#include <boost/filesystem/path.hpp>

#include "io/reference/reference_genome.hpp"
#include "basics/genomic_region.hpp"

namespace octopus { namespace io {

// Requires reference access to get contig sizes for partially specified regions (e.g. "4")
GenomicRegion parse_region(std::string region, const ReferenceGenome& reference);

enum class NonreferenceContigPolicy { exception, ignore };

std::deque<GenomicRegion> extract_regions(const boost::filesystem::path& file_path,
                                          const ReferenceGenome& reference,
                                          NonreferenceContigPolicy policy = NonreferenceContigPolicy::exception);

} // namespace io
} // namespace octopus

#endif
