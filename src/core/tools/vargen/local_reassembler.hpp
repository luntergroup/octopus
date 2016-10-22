// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef local_reassembler_hpp
#define local_reassembler_hpp

#include <vector>
#include <cstddef>
#include <functional>
#include <memory>

#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "concepts/mappable.hpp"
#include "core/types/variant.hpp"
#include "variant_generator.hpp"
#include "utils/assembler.hpp"

namespace octopus {

class ReferenceGenome;

namespace coretools {

class LocalReassembler : public VariantGenerator
{
public:
    struct Options
    {
        std::vector<unsigned> kmer_sizes              = {10, 25, 35};
        unsigned num_fallbacks                        = 6;
        unsigned fallback_interval_size               = 10;
        GenomicRegion::Size bin_size                  = 1000;
        GenomicRegion::Size bin_overlap               = 0;
        AlignedRead::BaseQuality mask_threshold       = 0;
        unsigned min_hard_prune_weight                = 1;
        double min_mean_path_weight                   = 2.0;
        unsigned max_paths                            = 10;
        Variant::MappingDomain::Size max_variant_size = 5000;
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
    
    GenomicRegion::Size bin_size_;
        
    GenomicRegion::Size bin_size_, bin_overlap_;
    std::deque<Bin> bins_;
    std::deque<NucleotideSequence> masked_sequence_buffer_;
    
    AlignedRead::BaseQuality mask_threshold_;
    unsigned min_hard_prune_weight_;
    double min_mean_path_weight_;
    unsigned max_paths_;
    Variant::MappingDomain::Size max_variant_size_;
    
    void prepare_bins_to_insert(const AlignedRead& read);
    unsigned try_assemble_with_defaults(const Bin& bin, std::deque<Variant>& result);
    void try_assemble_with_fallbacks(const Bin& bin, std::deque<Variant>& result);
    GenomicRegion propose_assembler_region(const GenomicRegion& input_region, unsigned kmer_size) const;
    bool assemble_bin(unsigned kmer_size, const Bin& bin, std::deque<Variant>& result) const;
    bool try_assemble_region(Assembler& assembler,
                             const NucleotideSequence& reference_sequence,
                             const GenomicRegion& reference_region,
                             std::deque<Variant>& result) const;
};

} // namespace coretools
} // namespace octopus

#endif
