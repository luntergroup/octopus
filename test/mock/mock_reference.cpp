// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mock_reference.hpp"

#include <iterator>
#include <algorithm>
#include <stdexcept>

namespace octopus { namespace test {
    
std::unique_ptr<ReferenceReader> MockReference::do_clone() const
{
    return std::make_unique<MockReference>(*this);
}

bool MockReference::do_is_open() const noexcept
{
    return true;
}

std::string MockReference::do_fetch_reference_name() const
{
    return "mock";
}

std::vector<MockReference::ContigName> MockReference::do_fetch_contig_names() const
{
    std::vector<ContigName> result {};
    result.reserve(mock_contigs_.size());
    std::transform(std::cbegin(mock_contigs_), std::cend(mock_contigs_), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    return result;
}

MockReference::GenomicSize MockReference::do_fetch_contig_size(const ContigName& contig) const
{
    return static_cast<GenomicSize>(mock_contigs_.at(contig).size());
}

MockReference::GeneticSequence MockReference::do_fetch_sequence(const GenomicRegion& region) const
{
    const auto& contig = mock_contigs_.at(contig_name(region));
    
    auto first = std::next(std::cbegin(contig), region.begin());
    
    if (region.end() > contig.size()) {
        throw std::runtime_error {"MockReference: out of bounds"};
    }
    
    auto last  = std::next(std::cbegin(contig), region.end());
    
    return {first, last};
}

ReferenceGenome make_mock_reference()
{
    return ReferenceGenome {std::make_unique<MockReference>()};
}

} // namespace test
} // namespace octopus
