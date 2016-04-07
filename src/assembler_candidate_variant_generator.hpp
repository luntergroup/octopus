//
//  assembler_candidate_variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__assembler_candidate_variant_generator__
#define __Octopus__assembler_candidate_variant_generator__

#include <vector>
#include <cstddef>
#include <functional>

#include "i_candidate_variant_generator.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
#include "assembler.hpp"

class ReferenceGenome;
class AlignedRead;
class GenomicRegion;
class Variant;

namespace Octopus
{
class AssemblerCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    using QualityType = AlignedRead::QualityType;
    using SizeType    = GenomicRegion::SizeType;
    
    AssemblerCandidateVariantGenerator() = delete;
    explicit AssemblerCandidateVariantGenerator(const ReferenceGenome& reference,
                                                std::vector<unsigned> kmer_sizes,
                                                QualityType min_base_quality = 10,
                                                unsigned min_supporting_reads = 2,
                                                SizeType max_variant_size = 500);
    ~AssemblerCandidateVariantGenerator() override = default;
    
    AssemblerCandidateVariantGenerator(const AssemblerCandidateVariantGenerator&)            = default;
    AssemblerCandidateVariantGenerator& operator=(const AssemblerCandidateVariantGenerator&) = default;
    AssemblerCandidateVariantGenerator(AssemblerCandidateVariantGenerator&&)                 = default;
    AssemblerCandidateVariantGenerator& operator=(AssemblerCandidateVariantGenerator&&)      = default;
    
    bool requires_reads() const noexcept override;
    void add_read(const AlignedRead& read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first,
                   std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                   MappableFlatMultiSet<AlignedRead>::const_iterator last) override;
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
    void clear() override;
    
private:
    using SequenceType = AlignedRead::SequenceType;
    
    struct Bin : public Mappable<Bin>
    {
        explicit Bin(GenomicRegion region);
        
        const GenomicRegion& get_region() const noexcept;
        
        void insert(const AlignedRead& read);
        void insert(const SequenceType& sequence);
        
        void clear() noexcept;
        
        GenomicRegion region;
        std::deque<std::reference_wrapper<const SequenceType>> read_sequences;
    };
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    std::vector<unsigned> default_kmer_sizes_;
    std::vector<unsigned> fallback_kmer_sizes_;
    
    SizeType bin_size_;
    std::deque<Bin> bins_;
    std::deque<SequenceType> masked_sequence_buffer_;
    
    QualityType min_base_quality_;
    unsigned min_supporting_reads_;
    SizeType max_variant_size_;
    
    void prepare_bins_to_insert(const AlignedRead& read);
    
    GenomicRegion propose_assembler_region(const GenomicRegion& input_region, unsigned kmer_size) const;
    bool assemble_bin(unsigned kmer_size, const Bin& bin, std::vector<Variant>& result) const;
    bool try_assemble_region(Assembler& assembler,
                             const SequenceType& reference_sequence,
                             const GenomicRegion& reference_region,
                             std::vector<Variant>& result) const;
};

} // namespace Octopus

#endif /* defined(__Octopus__assembler_candidate_variant_generator__) */
