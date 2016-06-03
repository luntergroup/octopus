//
//  candidate_generator_builder.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "candidate_generator_builder.hpp"

#include <algorithm>
#include <iterator>
#include <utility>

#include "common.hpp"
#include "logging.hpp"

namespace Octopus
{
    CandidateGeneratorBuilder::CandidateGeneratorBuilder()
    :
    generators_        {},
    parameters_        {},
    generator_factory_ {generate_factory()}
    {}
    
    CandidateGeneratorBuilder::CandidateGeneratorBuilder(const CandidateGeneratorBuilder& other)
    :
    generators_        {other.generators_},
    parameters_        {other.parameters_},
    generator_factory_ {generate_factory()}
    {}
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::operator=(const CandidateGeneratorBuilder& other)
    {
        generators_        = other.generators_;
        parameters_        = other.parameters_;
        generator_factory_ = generate_factory();
        return *this;
    }
    
    CandidateGeneratorBuilder::CandidateGeneratorBuilder(CandidateGeneratorBuilder&& other)
    :
    generators_        {std::move(other.generators_)},
    parameters_        {std::move(other.parameters_)},
    generator_factory_ {generate_factory()}
    {}
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::operator=(CandidateGeneratorBuilder&& other)
    {
        using std::swap;
        swap(generators_, other.generators_);
        swap(parameters_, other.parameters_);
        generator_factory_ = generate_factory();
        return *this;
    }
    
    unsigned CandidateGeneratorBuilder::num_generators() const noexcept
    {
        return static_cast<unsigned>(generators_.size());
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::add_generator(const Generator type)
    {
        if (std::find(std::cbegin(generators_), std::cend(generators_), type) == std::cend(generators_)) {
            generators_.push_back(type);
        }
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_reference(const ReferenceGenome& reference)
    {
        parameters_.reference = reference;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_min_base_quality(const QualityType quality)
    {
        parameters_.min_base_quality = quality;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_min_supporting_reads(const unsigned num_reads)
    {
        parameters_.min_supporting_reads = num_reads;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_max_variant_size(const SizeType size)
    {
        parameters_.max_variant_size = size;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::add_kmer_size(const unsigned kmer_size)
    {
        parameters_.kmer_sizes.push_back(kmer_size);
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_assembler_min_base_quality(const QualityType quality)
    {
        parameters_.min_assembler_base_quality = quality;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_variant_source(boost::filesystem::path variant_source)
    {
        parameters_.variant_source = std::make_shared<VcfReader>(std::move(variant_source));
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_variant_source(const std::shared_ptr<const VcfReader>& variant_source)
    {
        parameters_.variant_source = variant_source;
        return *this;
    }
    
    CandidateVariantGenerator CandidateGeneratorBuilder::build() const
    {
        CandidateVariantGenerator result {};
        
        for (const auto type : generators_) {
            result.register_generator(generator_factory_.at(type)());
        }
        
        return result;
    }
    
    // private methods
    
    CandidateGeneratorBuilder::GeneratorFactoryMap CandidateGeneratorBuilder::generate_factory() const
    {
        return GeneratorFactoryMap {
            {Generator::Alignment, [this] () {
                return std::make_unique<AlignmentCandidateVariantGenerator>(*parameters_.reference,
                                                                            parameters_.min_base_quality,
                                                                            parameters_.min_supporting_reads,
                                                                            parameters_.max_variant_size);
            }},
            {Generator::Assembler, [this] () {
                const auto quality = (parameters_.min_assembler_base_quality)
                        ? *parameters_.min_assembler_base_quality : parameters_.min_base_quality;
                return std::make_unique<AssemblerCandidateVariantGenerator>(*parameters_.reference,
                                                                            parameters_.kmer_sizes,
                                                                            quality,
                                                                            parameters_.min_supporting_reads,
                                                                            parameters_.max_variant_size);
            }},
            {Generator::External, [this] () {
                return std::make_unique<ExternalCandidateVariantGenerator>(parameters_.variant_source);
            }},
            {Generator::Online, [this] () {
                return std::make_unique<OnlineCandidateVariantGenerator>(*parameters_.reference,
                                                                         parameters_.max_variant_size);
            }},
            {Generator::Random, [this] () {
                return std::make_unique<RandomCandidateVariantGenerator>(*parameters_.reference);
            }}
        };
    }
} // namespace Octopus
