//
//  read_reader.h
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_reader__
#define __Octopus__read_reader__

#include <cstddef>
#include <iterator>
#include <memory>

#include <set>

#include "read_reader_implementor.h"
#include "genomic_region.h"
#include "aligned_read.h"

class ReadReader
{
public:
    ReadReader() = delete;
    ReadReader(const std::string& read_file_path);
    
    ReadReader(const ReadReader&)            = delete;
    ReadReader& operator=(const ReadReader&) = delete;
    ReadReader(ReadReader&&)                 = default;
    ReadReader& operator=(ReadReader&&)      = default;
    
    std::set<AlignedRead> fetch_reads(const GenomicRegion& a_region);
    
private:
    std::unique_ptr<IReadReaderImplementor> the_implementation_;
};

inline std::set<AlignedRead> ReadReader::fetch_reads(const GenomicRegion& a_region)
{
    return the_implementation_->fetch_reads(a_region);
}

#endif /* defined(__Octopus__read_reader__) */
