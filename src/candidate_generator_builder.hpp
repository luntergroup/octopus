//
//  candidate_generator_builder.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef candidate_generator_builder_hpp
#define candidate_generator_builder_hpp

#include <deque>
#include <functional>
#include <memory>
#include <unordered_map>
#include <cstddef>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "candidate_generators.hpp"
#include "reference_genome.hpp"
#include "vcf_reader.hpp"

namespace Octopus
{
    class CandidateGeneratorBuilder
    {
    public:
        using QualityType = AlignedRead::QualityType;
        using SizeType    = GenomicRegion::SizeType;
        
        enum class Generator { Alignment, Assembler, External, Online, Random };
        
        CandidateGeneratorBuilder();
        ~CandidateGeneratorBuilder() = default;
        
        CandidateGeneratorBuilder(const CandidateGeneratorBuilder&);
        CandidateGeneratorBuilder& operator=(const CandidateGeneratorBuilder&);
        CandidateGeneratorBuilder(CandidateGeneratorBuilder&&);
        CandidateGeneratorBuilder& operator=(CandidateGeneratorBuilder&&);
        
        unsigned num_generators() const noexcept;
        
        CandidateGeneratorBuilder& add_generator(Generator type);
        CandidateGeneratorBuilder& set_reference(const ReferenceGenome& reference);
        CandidateGeneratorBuilder& set_min_base_quality(QualityType quality);
        CandidateGeneratorBuilder& set_min_supporting_reads(unsigned num_reads);
        CandidateGeneratorBuilder& set_max_variant_size(SizeType size);
        
        CandidateGeneratorBuilder& add_kmer_size(unsigned kmer_size);
        CandidateGeneratorBuilder& set_assembler_min_base_quality(QualityType quality);
        
        CandidateGeneratorBuilder& set_variant_source(boost::filesystem::path variant_source);
        CandidateGeneratorBuilder& set_variant_source(const std::shared_ptr<VcfReader>& variant_source);
        
        CandidateVariantGenerator build() const;
        
    private:
        struct Parameters
        {
            Parameters()  = default;
            ~Parameters() = default;
            
            Parameters(const Parameters&)            = default;
            Parameters& operator=(const Parameters&) = default;
            Parameters(Parameters&&)                 = default;
            Parameters& operator=(Parameters&&)      = default;
            
            // common
            boost::optional<std::reference_wrapper<const ReferenceGenome>> reference;
            
            QualityType min_base_quality = 10;
            unsigned min_supporting_reads = 1;
            SizeType max_variant_size = 500;
            
            // assembler
            
            std::vector<unsigned> kmer_sizes;
            boost::optional<QualityType> min_assembler_base_quality;
            
            // external
            
            std::shared_ptr<VcfReader> variant_source;
            
            // online
            
            // random
        };
        
        struct GeneratorTypeHash
        {
            std::size_t operator()(Generator type) const { return static_cast<std::size_t>(type); }
        };
        
        using GeneratorFactoryMap = std::unordered_map<Generator,
                std::function<std::unique_ptr<ICandidateVariantGenerator>()>, GeneratorTypeHash>;
        
        std::deque<Generator> generators_;
        Parameters parameters_;
        GeneratorFactoryMap generator_factory_;
        
        GeneratorFactoryMap generate_factory() const;
    };
} // namespace Octopus

#endif /* candidate_generator_builder_hpp */
