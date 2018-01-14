// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "active_region_generator.hpp"

namespace octopus { namespace coretools {

void ActiveRegionGenerator::add_read(const SampleName& sample, const AlignedRead& read)
{

}

std::vector<GenomicRegion> ActiveRegionGenerator::generate(const GenomicRegion& region, const std::string& generator) const
{
    return {};
}

void ActiveRegionGenerator::clear() noexcept
{

}

} // namespace coretools
} // namespace octopus
