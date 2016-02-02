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

#include <iostream> // DEBUG

namespace Octopus
{
    CandidateGeneratorBuilder::CandidateGeneratorBuilder()
    :
    generators_           {},
    reference_            {},
    min_snp_base_quality_ {},
    min_supporting_reads_ {},
    max_variant_size_     {},
    kmer_size_            {},
    variant_source_       {},
    generator_factory_    {generate_factory()}
    {}
    
    CandidateGeneratorBuilder::CandidateGeneratorBuilder(const CandidateGeneratorBuilder& other)
    :
    generators_           {other.generators_},
    reference_            {other.reference_},
    min_snp_base_quality_ {other.min_snp_base_quality_},
    min_supporting_reads_ {other.min_supporting_reads_},
    max_variant_size_     {other.max_variant_size_},
    kmer_size_            {other.kmer_size_},
    variant_source_       {other.variant_source_},
    generator_factory_    {generate_factory()}
    {}
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::operator=(const CandidateGeneratorBuilder& other)
    {
        generators_           = other.generators_;
        reference_            = other.reference_;
        min_snp_base_quality_ = other.min_snp_base_quality_;
        min_supporting_reads_ = other.min_supporting_reads_;
        max_variant_size_     = other.max_variant_size_;
        kmer_size_            = other.kmer_size_;
        variant_source_       = other.variant_source_;
        generator_factory_    = generate_factory();
        return *this;
    }
    
    CandidateGeneratorBuilder::CandidateGeneratorBuilder(CandidateGeneratorBuilder&& other)
    :
    generators_           {std::move(other.generators_)},
    reference_            {std::move(other.reference_)},
    min_snp_base_quality_ {std::move(other.min_snp_base_quality_)},
    min_supporting_reads_ {std::move(other.min_supporting_reads_)},
    max_variant_size_     {std::move(other.max_variant_size_)},
    kmer_size_            {std::move(other.kmer_size_)},
    variant_source_       {std::move(other.variant_source_)},
    generator_factory_    {generate_factory()}
    {}
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::operator=(CandidateGeneratorBuilder&& other)
    {
        using std::swap;
        swap(generators_,           other.generators_);
        swap(reference_,            other.reference_);
        swap(min_snp_base_quality_, other.min_snp_base_quality_);
        swap(min_supporting_reads_, other.min_supporting_reads_);
        swap(max_variant_size_,     other.max_variant_size_);
        swap(kmer_size_,            other.kmer_size_);
        swap(variant_source_,       other.variant_source_);
        
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
        reference_ = reference;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_min_snp_base_quality(const QualityType quality)
    {
        min_snp_base_quality_ = quality;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_min_supporting_reads(const unsigned num_reads)
    {
        min_supporting_reads_ = num_reads;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_max_variant_size(const SizeType size)
    {
        max_variant_size_ = size;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_kmer_size(const unsigned kmer_size)
    {
        kmer_size_ = kmer_size;
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_variant_source(boost::filesystem::path variant_source)
    {
        variant_source_ = std::make_shared<VcfReader>(std::move(variant_source));
        return *this;
    }
    
    CandidateGeneratorBuilder& CandidateGeneratorBuilder::set_variant_source(const std::shared_ptr<VcfReader>& variant_source)
    {
        variant_source_ = variant_source;
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
                return std::make_unique<AlignmentCandidateVariantGenerator>(*reference_,
                                                                            min_snp_base_quality_,
                                                                            min_supporting_reads_,
                                                                            max_variant_size_);
            }},
            {Generator::Assembler, [this] () {
                return std::make_unique<AssemblerCandidateVariantGenerator>(*reference_,
                                                                            *kmer_size_,
                                                                            max_variant_size_);
            }},
            {Generator::External, [this] () {
                return std::make_unique<ExternalCandidateVariantGenerator>(variant_source_);
            }},
            {Generator::Online, [this] () {
                return std::make_unique<OnlineCandidateVariantGenerator>(*reference_,
                                                                         max_variant_size_);
            }},
            {Generator::Random, [this] () {
                return std::make_unique<RandomCandidateVariantGenerator>(*reference_);
            }}
        };
    }
} // namespace Octopus
