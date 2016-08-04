// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef Octopus_test_utils_h
#define Octopus_test_utils_h

#include <cmath>
#include <string>
#include <unordered_map>
#include <initializer_list>
#include <algorithm>
#include <iterator>

#include <core/types/haplotype.hpp>

namespace octopus {

class GenomicRegion;
class ReferenceGenome;

namespace test {

template <typename RealType>
bool is_close_to_one(RealType x)
{
    return std::abs(x - 1) < 0.0000000000001;
}

inline std::string get_label(const std::unordered_map<Haplotype, std::string>& labels, const Haplotype& haplotype)
{
    return (labels.count(haplotype) > 0) ? labels.at(haplotype) : "other";
}

inline void sort(std::vector<Haplotype>& haplotypes)
{
    std::sort(std::begin(haplotypes), std::end(haplotypes));
}

template <typename A>
Haplotype make_haplotype(const ReferenceGenome& reference, const GenomicRegion& region,
                         std::initializer_list<A> alleles)
{
    Haplotype::Builder hb {region, reference};
    
    for (const auto& allele : alleles) hb.push_back(allele);
    
    return hb.build();
}

} // namespace test
} // namespace octopus

#endif
