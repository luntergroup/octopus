// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "ploidies.hpp"

#include <utility>
#include <iterator>
#include <algorithm>

namespace octopus { namespace csr {

const std::string Ploidies::name_ {"Ploidies"};

Ploidies::Ploidies(const PloidyMap& ploidies, const GenomicRegion& region, const std::vector<SampleName>& samples)
{
    const auto local_ploidies = get_ploidies(samples, region.contig_name(), ploidies);
    ploidies_.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::cbegin(local_ploidies),
                   std::inserter(ploidies_, std::begin(ploidies_)),
                   [] (const auto& sample, auto ploidy) { return std::make_pair(sample, ploidy); });
}

Facet::ResultType Ploidies::do_get() const
{
    return std::cref(ploidies_);
}

} // namespace csr
} // namespace octopus
