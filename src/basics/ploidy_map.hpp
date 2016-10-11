// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef ploidy_map_hpp
#define ploidy_map_hpp

#include <unordered_map>

#include "config/common.hpp"

namespace octopus {

class PloidyMap
{
public:
    PloidyMap(unsigned organism = 2);
    
    PloidyMap(const PloidyMap&)            = default;
    PloidyMap& operator=(const PloidyMap&) = default;
    PloidyMap(PloidyMap&&)                 = default;
    PloidyMap& operator=(PloidyMap&&)      = default;
    
    ~PloidyMap() = default;
    
    void set(const SampleName& sample, const ContigName& contig, unsigned ploidy);
    void set(const ContigName& contig, unsigned ploidy);
    
    unsigned organism() const noexcept;
    unsigned operator()(const SampleName& sample, const ContigName& contig) const noexcept;
    
private:
    using ContigPloidyMap = std::unordered_map<ContigName, unsigned>;
    
    unsigned organism_;
    ContigPloidyMap contigs_;
    std::unordered_map<SampleName, ContigPloidyMap> sample_contigs_;
};

} // namespace octopus

#endif
