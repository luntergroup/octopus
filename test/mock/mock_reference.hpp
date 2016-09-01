// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mock_reference_hpp
#define mock_reference_hpp

#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include "io/reference/reference_reader.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus { namespace test { namespace mock {

using octopus::io::ReferenceReader;

class MockReference : public ReferenceReader
{
public:
    using ContigName      = ReferenceReader::ContigName;
    using GenomicSize     = ReferenceReader::GenomicSize;
    using GeneticSequence = ReferenceReader::GeneticSequence;
    
    MockReference();
    
    MockReference(const MockReference&)            = default;
    MockReference& operator=(const MockReference&) = default;
    MockReference(MockReference&&)                 = default;
    MockReference& operator=(MockReference&&)      = default;
    
private:
    std::unique_ptr<ReferenceReader> do_clone() const override;
    bool do_is_open() const noexcept override;
    std::string do_fetch_reference_name() const override;
    std::vector<ContigName> do_fetch_contig_names() const override;
    GenomicSize do_fetch_contig_size(const ContigName& contig) const override;
    GeneticSequence do_fetch_sequence(const GenomicRegion& region) const override;
    
    std::unordered_map<ContigName, GeneticSequence> mock_contigs_;
};
    
ReferenceGenome make_reference();
    
} // namespace mock
} // namespace test
} // namespace octopus

#endif
