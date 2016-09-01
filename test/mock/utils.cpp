// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "utils.hpp"

#include <utility>

#include "io/region/region_parser.hpp"

namespace octopus { namespace debug {

Allele make_allele(const std::string& region, std::string sequence, const ReferenceGenome& reference)
{
    return Allele {io::parse_region(region, reference), std::move(sequence)};
}

Allele make_reference_allele(const std::string& region, const ReferenceGenome& reference)
{
    return make_reference_allele(io::parse_region(region, reference), reference);
}

Variant make_variant(const std::string& region_str, Variant::NucleotideSequence alt_sequence,
                     const ReferenceGenome& reference)
{
    return make_variant(Allele {io::parse_region(region_str, reference), std::move(alt_sequence)}, reference);
}

Haplotype make_haplotype(const std::string& str, const GenomicRegion& region,
                         const ReferenceGenome& reference)
{
    if (str.size() < 3 || str.front() != '<' || str.back() != '>') {
        throw std::runtime_error {"make_haplotype: bad input"};
    }
    
    Haplotype::Builder hb {region, reference};
    
    if (str.size() == 3) {
        return hb.build(); // reference
    }
    
    std::size_t pos {3};
    
    while (pos < str.size()) {
        auto i = str.find(' ', pos);
        auto j = str.find('}', i + 1);
        hb.push_back(make_allele(str.substr(pos, i - pos), str.substr(i + 1, j - i - 1), reference));
        pos = j + 3;
    }
    
    return hb.build();
}

Haplotype make_haplotype(const std::string& str, const std::string& region,
                         const ReferenceGenome& reference)
{
    return make_haplotype(str, io::parse_region(region, reference), reference);
}

Genotype<Haplotype> make_genotype(const std::string& str, const GenomicRegion& region,
                                  const ReferenceGenome& reference)
{
    Genotype<Haplotype> result {};
    
    //std::size_t pos {0};
    
    //        while (pos != str.size()) {
    //
    //        }
    
    return result;
}

Genotype<Haplotype> make_genotype(const std::string& str, const std::string& region,
                                  const ReferenceGenome& reference)
{
    return make_genotype(str, io::parse_region(region, reference), reference);
}

} // namespace debug
} // namespace octopus
