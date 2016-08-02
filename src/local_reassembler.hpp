//
//  local_reassembler.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__local_reassembler__
#define __Octopus__local_reassembler__

#include <vector>
#include <cstddef>
#include <functional>
#include <memory>

#include "variant_generator.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
#include "assembler.hpp"
#include "variant.hpp"

class ReferenceGenome;
class AlignedRead;
class GenomicRegion;

namespace octopus { namespace coretools {

class LocalReassembler : public VariantGenerator
{
public:
    struct Options
    {
        std::vector<unsigned> kmer_sizes           = {10, 25, 35};
        AlignedRead::BaseQuality mask_threshold    = 0;
        unsigned min_supporting_reads              = 2;
        Variant::RegionType::Size max_variant_size = 500;
    };
    
    LocalReassembler() = delete;
    
    LocalReassembler(const ReferenceGenome& reference, Options options);
    
    LocalReassembler(const LocalReassembler&)            = default;
    LocalReassembler& operator=(const LocalReassembler&) = default;
    LocalReassembler(LocalReassembler&&)                 = default;
    LocalReassembler& operator=(LocalReassembler&&)      = default;
    
    ~LocalReassembler() override = default;
    
private:
    using VariantGenerator::VectorIterator;
    using VariantGenerator::FlatSetIterator;
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    
    bool do_requires_reads() const noexcept override;
    
    void do_add_read(const AlignedRead& read) override;
    void do_add_reads(VectorIterator first, VectorIterator last) override;
    void do_add_reads(FlatSetIterator first, FlatSetIterator last) override;
    
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    
    void do_clear() noexcept override;
    
    std::string name() const override;
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    
    struct Bin : public Mappable<Bin>
    {
        Bin(GenomicRegion region);
        
        const GenomicRegion& mapped_region() const noexcept;
        
        void insert(const AlignedRead& read);
        void insert(const NucleotideSequence& sequence);
        
        void clear() noexcept;
        bool empty() const noexcept;
        
        GenomicRegion region;
        std::deque<std::reference_wrapper<const NucleotideSequence>> read_sequences;
    };
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    std::vector<unsigned> default_kmer_sizes_;
    std::vector<unsigned> fallback_kmer_sizes_;
    
    ContigRegion::Size bin_size_;
    std::deque<Bin> bins_;
    std::deque<NucleotideSequence> masked_sequence_buffer_;
    
    AlignedRead::BaseQuality mask_threshold_;
    unsigned min_supporting_reads_;
    Variant::RegionType::Size max_variant_size_;
    
    void prepare_bins_to_insert(const AlignedRead& read);
    
    GenomicRegion propose_assembler_region(const GenomicRegion& input_region, unsigned kmer_size) const;
    bool assemble_bin(unsigned kmer_size, const Bin& bin, std::deque<Variant>& result) const;
    bool try_assemble_region(Assembler& assembler,
                             const NucleotideSequence& reference_sequence,
                             const GenomicRegion& reference_region,
                             std::deque<Variant>& result) const;
};

} // namespace coretools
} // namespace octopus

#endif /* defined(__Octopus__local_reassembler__) */
