// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reference_reader_hpp
#define reference_reader_hpp

#include <string>
#include <vector>
#include <cstdint>
#include <memory>

#include "basics/genomic_region.hpp"

namespace octopus { namespace io {

class ReferenceReader
{
public:
    using ContigName      = GenomicRegion::ContigName;
    using GenomicSize     = GenomicRegion::Size;
    using GeneticSequence = std::string;
    
    virtual ~ReferenceReader() = default;
    
    std::unique_ptr<ReferenceReader> clone() const
    {
        return do_clone();
    }
    
    auto is_open() const noexcept
    {
        return do_is_open();
    }
    
    auto fetch_reference_name() const
    {
        return do_fetch_reference_name();
    }
    
    auto fetch_contig_names() const
    {
        return do_fetch_contig_names();
    }
    
    auto fetch_contig_size(const ContigName& contig) const
    {
        return do_fetch_contig_size(contig);
    }
    
    GeneticSequence fetch_sequence(const GenomicRegion& region) const
    {
        return do_fetch_sequence(region);
    }
    
private:
    virtual std::unique_ptr<ReferenceReader> do_clone() const = 0;
    virtual bool do_is_open() const noexcept = 0;
    virtual std::string do_fetch_reference_name() const = 0;
    virtual std::vector<ContigName> do_fetch_contig_names() const = 0;
    virtual GenomicSize do_fetch_contig_size(const ContigName& contig) const = 0;
    virtual GeneticSequence do_fetch_sequence(const GenomicRegion& region) const = 0;
};

} // namespace io
} // namespace octopus

#endif
