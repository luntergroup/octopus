// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef memory_footprint_hpp
#define memory_footprint_hpp

#include <cstddef>
#include <string>
#include <functional>
#include <iosfwd>

#include <boost/optional.hpp>

#include "concepts/comparable.hpp"

namespace octopus {

class MemoryFootprint : public Comparable<MemoryFootprint>
{
public:
    MemoryFootprint() = default;
    
    constexpr MemoryFootprint(std::size_t num_bytes) noexcept : num_bytes_ {num_bytes} {}
    
    MemoryFootprint(const MemoryFootprint&)            = default;
    MemoryFootprint& operator=(const MemoryFootprint&) = default;
    MemoryFootprint(MemoryFootprint&&)                 = default;
    MemoryFootprint& operator=(MemoryFootprint&&)      = default;
    
    ~MemoryFootprint() = default;
    
    constexpr std::size_t num_bytes() const noexcept { return num_bytes_; }
    
private:
    std::size_t num_bytes_;
};

bool operator==(const MemoryFootprint& lhs, const MemoryFootprint& rhs) noexcept;
bool operator<(const MemoryFootprint& lhs, const MemoryFootprint& rhs) noexcept;

std::ostream& operator<<(std::ostream& os, MemoryFootprint footprint);
std::istream& operator>>(std::istream& is, MemoryFootprint& result);

boost::optional<MemoryFootprint> parse_footprint(std::string footprint_str);

} // namespace octopus

namespace std {

template <> struct hash<octopus::MemoryFootprint>
{
    size_t operator()(const octopus::MemoryFootprint& fp) const noexcept
    {
        return hash<decltype(fp.num_bytes())>()(fp.num_bytes());
    }
};

} // namespace std

#endif
