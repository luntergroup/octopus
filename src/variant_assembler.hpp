//
//  variant_assembler.hpp
//  Octopus
//
//  Created by Daniel Cooke on 24/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_assembler__
#define __Octopus__variant_assembler__

#include "kmer_graph.hpp"

#include <string>
#include <vector>

#include "storage_policies.hpp"

class Variant;
class AlignedRead;
class GenomicRegion;

class VariantAssembler
{
public:
    VariantAssembler() = delete;
    explicit VariantAssembler(unsigned k);
    ~VariantAssembler() = default;
    
    VariantAssembler(const VariantAssembler&)            = default;
    VariantAssembler& operator=(const VariantAssembler&) = default;
    VariantAssembler(VariantAssembler&&)                 = default;
    VariantAssembler& operator=(VariantAssembler&&)      = default;
    
    void add_read(const AlignedRead& read);
    void add_reference_sequence(const GenomicRegion& region, const std::string& sequence);
    std::vector<Variant> get_variants(const GenomicRegion& region);
    void clear() noexcept;
    
private:
    enum class Colour {Reference, Read};
    
    KmerGraph<Colour, policies::StoreStringReference> de_bruijn_graph_;
};

#endif /* defined(__Octopus__variant_assembler__) */
