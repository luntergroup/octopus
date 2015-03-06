//
//  read_filter.h
//  Octopus
//
//  Created by Daniel Cooke on 06/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filter__
#define __Octopus__read_filter__

#include <functional>

class AlignedRead;

class ReadFilter
{
public:
    ReadFilter() = default;
    ~ReadFilter() = default;
    
    ReadFilter(const ReadFilter&)            = default;
    ReadFilter& operator=(const ReadFilter&) = default;
    ReadFilter(ReadFilter&&)                 = default;
    ReadFilter& operator=(ReadFilter&&)      = default;
    
    void register_filter(std::function<AlignedRead(AlignedRead)>);
    unsigned num_filters() const noexcept;
    AlignedRead filter_read(const AlignedRead& a_read);
    
private:
    
};

#endif /* defined(__Octopus__read_filter__) */
