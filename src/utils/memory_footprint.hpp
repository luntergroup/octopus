// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef memory_footprint_hpp
#define memory_footprint_hpp

#include <cstddef>
#include <string>
#include <iosfwd>

#include <boost/optional.hpp>

namespace octopus {

class MemoryFootprint
{
public:
    MemoryFootprint() = default;
    
    MemoryFootprint(std::size_t num_bytes) noexcept;
    
    MemoryFootprint(const MemoryFootprint&)            = default;
    MemoryFootprint& operator=(const MemoryFootprint&) = default;
    MemoryFootprint(MemoryFootprint&&)                 = default;
    MemoryFootprint& operator=(MemoryFootprint&&)      = default;
    
    ~MemoryFootprint() = default;
        
    std::size_t num_bytes() const noexcept;
    
private:
    std::size_t num_bytes_;
};

std::ostream& operator<<(std::ostream& os, MemoryFootprint footprint);
std::istream& operator>>(std::istream& is, MemoryFootprint& result);

boost::optional<MemoryFootprint> parse_footprint(std::string str);

} // namespace octopus

#endif
