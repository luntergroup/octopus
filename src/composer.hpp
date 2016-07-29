//
//  composer.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__composer__
#define __Octopus__composer__

#include <vector>
#include <memory>
#include <cstddef>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "variant.hpp"
#include "candidate_variant_generator.hpp"
#include "reference_genome.hpp"
#include "vcf_reader.hpp"

class GenomicRegion;
class AlignedRead;

namespace octopus { namespace core { namespace generators
{
class Composer : public CandidateVariantGenerator
{
public:
    class Builder;
    
    Composer() = default;
    
    Composer(const Composer&)            = delete;
    Composer& operator=(const Composer&) = delete;
    Composer(Composer&&)                 = default;
    Composer& operator=(Composer&&)      = default;
    
    ~Composer() override = default;
    
    void register_generator(std::unique_ptr<CandidateVariantGenerator> generator);
    
    bool requires_reads() const noexcept override;
    
    void add_read(const AlignedRead& read) override;
    void add_reads(std::vector<AlignedRead>::const_iterator first,
                   std::vector<AlignedRead>::const_iterator last) override;
    void add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                   MappableFlatMultiSet<AlignedRead>::const_iterator last) override;
    
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
    
    void reserve(std::size_t n) override;
    void clear() override;
    
private:
    std::vector<std::unique_ptr<CandidateVariantGenerator>> generators_;
};

class Composer::Builder
{
public:
    using BaseQuality = AlignedRead::BaseQuality;
    using Position    = GenomicRegion::Position;
    
    enum class Generator { Alignment, Assembler, External, Online, Random };
    
    Builder();
    
    Builder(const Builder&);
    Builder& operator=(const Builder&);
    Builder(Builder&&);
    Builder& operator=(Builder&&);
    
    ~Builder() = default;
    
    unsigned num_generators() const noexcept;
    
    Builder& add_generator(Generator type);
    Builder& set_reference(const ReferenceGenome& reference);
    Builder& set_min_base_quality(BaseQuality min);
    Builder& set_min_supporting_reads(unsigned num_reads);
    Builder& set_max_variant_size(Position size);
    
    Builder& add_kmer_size(unsigned kmer_size);
    Builder& set_assembler_min_base_quality(BaseQuality min);
    
    Builder& set_variant_source(boost::filesystem::path variant_source);
    Builder& set_variant_source(const std::shared_ptr<const VcfReader>& variant_source);
    
    Composer build() const;
    
private:
    struct Parameters
    {
        // common
        boost::optional<std::reference_wrapper<const ReferenceGenome>> reference;
        
        BaseQuality min_base_quality = 10;
        unsigned min_supporting_reads = 1;
        Position max_variant_size = 500;
        
        // assembler
        std::vector<unsigned> kmer_sizes;
        boost::optional<BaseQuality> min_assembler_base_quality;
        
        // external
        std::shared_ptr<const VcfReader> variant_source;
    };
    
    struct GeneratorTypeHash
    {
        std::size_t operator()(Generator type) const { return static_cast<std::size_t>(type); }
    };
    
    using GeneratorFactoryMap = std::unordered_map<Generator,
    std::function<std::unique_ptr<CandidateVariantGenerator>()>, GeneratorTypeHash>;
    
    std::deque<Generator> generators_;
    
    Parameters parameters_;
    
    GeneratorFactoryMap generator_factory_;
    
    GeneratorFactoryMap generate_factory() const;
};
} // namespace generators
} // namespace core
} // namespace octopus

#endif /* defined(__Octopus__variant_candidate_generator__) */
