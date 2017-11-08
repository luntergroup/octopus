// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef utils_hpp
#define utils_hpp

#include <string>

#include "io/reference/reference_genome.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"

namespace octopus { namespace debug {
    
Allele make_allele(const std::string& region, std::string sequence, const ReferenceGenome& reference);

Allele make_reference_allele(const std::string& region, const ReferenceGenome& reference);

Variant make_variant(const std::string& region_str, Variant::NucleotideSequence alt_sequence,
                     const ReferenceGenome& reference);

Haplotype make_haplotype(const std::string& str, const GenomicRegion& region,
                         const ReferenceGenome& reference);

Haplotype make_haplotype(const std::string& str, const std::string& region,
                         const ReferenceGenome& reference);

Genotype<Haplotype> make_genotype(const std::string& str, const GenomicRegion& region,
                                  const ReferenceGenome& reference);

Genotype<Haplotype> make_genotype(const std::string& str, const std::string& region,
                                  const ReferenceGenome& reference);

} // namespace debug
} // namespace octopus

#endif
