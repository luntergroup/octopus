// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef ploidy_map_hpp
#define ploidy_map_hpp

#include <unordered_map>
#include <vector>

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
    
    bool is_autosome(const ContigName& contig) const noexcept;
    
    unsigned of(const SampleName& sample, const ContigName& contig) const noexcept;
    
private:
    using ContigPloidyMap = std::unordered_map<ContigName, unsigned>;
    
    unsigned organism_;
    ContigPloidyMap contigs_;
    std::unordered_map<ContigName, std::unordered_map<SampleName, unsigned>> allosomes_;
};

std::vector<unsigned> get_ploidies(const std::vector<SampleName>& samples, const ContigName& contig, const PloidyMap& ploidies);

unsigned get_min_ploidy(const std::vector<SampleName>& samples, const std::vector<ContigName>& contigs, const PloidyMap& ploidies);
unsigned get_max_ploidy(const std::vector<SampleName>& samples, const std::vector<ContigName>& contigs, const PloidyMap& ploidies);

} // namespace octopus

#endif
