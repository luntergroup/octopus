// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef local_reassembler_hpp
#define local_reassembler_hpp

#include <vector>
#include <cstddef>
#include <functional>
#include <memory>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "basics/mappable_reference_wrapper.hpp"
#include "concepts/mappable.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "core/types/variant.hpp"
#include "variant_generator.hpp"
#include "utils/assembler.hpp"

#include "utils/assembler_active_region_generator.hpp"

namespace octopus {

class ReferenceGenome;

namespace coretools {

class LocalReassembler : public VariantGenerator
{
public:
    struct Options
    {
        ExecutionPolicy execution_policy              = ExecutionPolicy::seq;
        std::vector<unsigned> kmer_sizes              = {10, 25, 35};
        unsigned num_fallbacks                        = 6;
        unsigned fallback_interval_size               = 10;
        GenomicRegion::Size bin_size                  = 1000;
        GenomicRegion::Size bin_overlap               = 0;
        AlignedRead::BaseQuality mask_threshold       = 0;
        unsigned min_kmer_observations                = 1;
        unsigned max_bubbles                          = 10;
        double min_bubble_score                       = 2.0;
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
    
    void do_add_read(const SampleName& sample, const AlignedRead& read) override;
    void do_add_reads(const SampleName& sample, VectorIterator first, VectorIterator last) override;
    void do_add_reads(const SampleName& sample, FlatSetIterator first, FlatSetIterator last) override;
    
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    
    void do_clear() noexcept override;
    
    std::string name() const override;
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    
    struct Bin : public Mappable<Bin>
    {
        Bin(GenomicRegion region);
        
        const GenomicRegion& mapped_region() const noexcept;
        
        void add(const AlignedRead& read);
        void add(const GenomicRegion& read_region, const NucleotideSequence& read_sequence);
        
        void clear() noexcept;
        bool empty() const noexcept;
        
        GenomicRegion region;
        boost::optional<ContigRegion> read_region;
        std::deque<std::reference_wrapper<const NucleotideSequence>> read_sequences;
    };
    
    ExecutionPolicy execution_policy_;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    std::vector<unsigned> default_kmer_sizes_, fallback_kmer_sizes_;
    
    MappableFlatMultiSet<MappableReferenceWrapper<const AlignedRead>> read_buffer_;
    
    GenomicRegion::Size bin_size_, bin_overlap_;
    std::deque<Bin> bins_;
    std::deque<NucleotideSequence> masked_sequence_buffer_;
    
    AlignedRead::BaseQuality mask_threshold_;
    unsigned min_kmer_observations_;
    unsigned max_bubbles_;
    double min_bubble_score_;
    Variant::MappingDomain::Size max_variant_size_;
    
    AssemblerActiveRegionGenerator active_region_generator_;
    
    void prepare_bins_to_insert(const AlignedRead& read);
    bool should_assemble_bin(const Bin& bin) const;
    void finalise_bins();
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
