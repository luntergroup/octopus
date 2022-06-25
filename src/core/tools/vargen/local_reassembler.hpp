// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef local_reassembler_hpp
#define local_reassembler_hpp

#include <vector>
#include <map>
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

namespace octopus {

class ReferenceGenome;

namespace coretools {

class LocalReassembler : public VariantGenerator
{
public:
    using ReadBaseCountMap = std::unordered_map<SampleName, unsigned>;
    using BubbleScoreSetter = std::function<double(const GenomicRegion&, const ReadBaseCountMap&)>;
    
    struct Options
    {
        enum class CyclicGraphTolerance { high, low, none };
        ExecutionPolicy execution_policy              = ExecutionPolicy::seq;
        std::vector<unsigned> kmer_sizes              = {10, 25, 35};
        unsigned num_fallbacks                        = 6;
        unsigned fallback_interval_size               = 10;
        GenomicRegion::Size bin_size                  = 1000;
        GenomicRegion::Size bin_overlap               = 0;
        AlignedRead::BaseQuality mask_threshold       = 0;
        unsigned min_kmer_observations                = 1;
        unsigned max_bubbles                          = 10;
        BubbleScoreSetter min_bubble_score            = [] (const GenomicRegion&, const ReadBaseCountMap&) { return 2.0; };
        Variant::MappingDomain::Size max_variant_size = 5000;
        CyclicGraphTolerance cycle_tolerance          = CyclicGraphTolerance::high;
        bool ignore_strand_bias                       = false;
    };
    
    LocalReassembler() = delete;
    
    LocalReassembler(const ReferenceGenome& reference, Options options);
    
    LocalReassembler(const LocalReassembler&)            = default;
    LocalReassembler& operator=(const LocalReassembler&) = default;
    LocalReassembler(LocalReassembler&&)                 = default;
    LocalReassembler& operator=(LocalReassembler&&)      = default;
    
    ~LocalReassembler() override = default;

private:
    using VariantGenerator::ReadVectorIterator;
    using VariantGenerator::ReadFlatSetIterator;
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    
    bool do_requires_reads() const noexcept override;
    
    void do_add_read(const SampleName& sample, const AlignedRead& read) override;
    void do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last) override;
    void do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last) override;
    
    std::vector<Variant> do_generate(const RegionSet& regions, OptionalThreadPool workers) const override;
    
    void do_clear() noexcept override;
    
    std::string name() const override;
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using SequenceBuffer     = std::deque<NucleotideSequence>;
    using ReadReference      = MappableReferenceWrapper<const AlignedRead>;
    using ReadBuffer         = MappableFlatMultiSet<ReadReference>;
    using ReadBufferMap      = std::map<SampleName, ReadBuffer>;
    
    struct Bin : public Mappable<Bin>
    {
        struct ReadData
        {
            std::reference_wrapper<const NucleotideSequence> sequence;
            std::reference_wrapper<const AlignedRead::BaseQualityVector> base_qualities;
            std::size_t sample_index;
        };
        using ReadDataStash = std::deque<ReadData>;
        
        Bin(GenomicRegion region);
        
        const GenomicRegion& mapped_region() const noexcept;
        
        void add(const AlignedRead& read, std::size_t sample_index);
        void add(const AlignedRead& read, const NucleotideSequence& masked_sequence, std::size_t sample_index);
        
        void clear() noexcept;
        std::size_t size() const noexcept;
        bool empty() const noexcept;
        
        GenomicRegion region;
        boost::optional<ContigRegion> read_region;
        ReadDataStash forward_read_sequences, reverse_read_sequences;
    };
    
    using BinList = std::deque<Bin>;
    
    enum class AssemblerStatus { success, partial_success, failed };
    
    ExecutionPolicy execution_policy_;
    std::reference_wrapper<const ReferenceGenome> reference_;
    std::vector<unsigned> default_kmer_sizes_, fallback_kmer_sizes_;
    mutable ReadBufferMap read_buffer_;
    GenomicRegion::Size max_bin_size_, max_bin_overlap_;
    AlignedRead::BaseQuality mask_threshold_;
    unsigned min_kmer_observations_;
    unsigned max_bubbles_;
    BubbleScoreSetter min_bubble_score_;
    Variant::MappingDomain::Size max_variant_size_;
    Options::CyclicGraphTolerance cycle_tolerance_;
    bool ignore_strand_bias_;
    
    void prepare_bins(const GenomicRegion& active_region, BinList& bins) const;
    bool should_assemble_bin(const Bin& bin) const;
    void finalise_bins(BinList& bins, const RegionSet& active_regions) const;
    unsigned try_assemble_with_defaults(const Bin& bin, std::deque<Variant>& result) const;
    void try_assemble_with_fallbacks(const Bin& bin, std::deque<Variant>& result) const;
    GenomicRegion propose_assembler_region(const GenomicRegion& input_region, unsigned kmer_size) const;
    void load(const Bin& bin, Assembler& assembler) const;
    void load(const Bin& bin, std::size_t sample_idx, Assembler& assembler) const;
    AssemblerStatus assemble_bin(unsigned kmer_size, const Bin& bin, std::deque<Variant>& result) const;
    AssemblerStatus try_assemble_region(Assembler& assembler, const NucleotideSequence& reference_sequence,
                                        const GenomicRegion& reference_region, std::deque<Variant>& result) const;
    double calculate_min_bubble_score(const GenomicRegion& bubble_region) const;
};

struct ConstantBubbleScoreSetter
{
    double operator()(const GenomicRegion&, const LocalReassembler::ReadBaseCountMap&) const noexcept { return min_score_; }
    ConstantBubbleScoreSetter(double min_score) : min_score_ {min_score} {}
private:
    double min_score_;
};

struct DepthBasedBubbleScoreSetter
{
    double operator()(const GenomicRegion& region, const LocalReassembler::ReadBaseCountMap& read_counts) const;
    DepthBasedBubbleScoreSetter(double min_score, double min_allele_frequency);
private:
    
    double min_score_, min_allele_frequency_;
};

} // namespace coretools
} // namespace octopus

#endif
