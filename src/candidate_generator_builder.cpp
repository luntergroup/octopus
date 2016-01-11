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
    generator_factory_ {
        {Generator::Alignment, [&] () {
            return std::make_unique<AlignmentCandidateVariantGenerator>(reference_.get(),
                                                                        min_snp_base_quality_,
                                                                        min_supporting_reads_,
                                                                        max_variant_size_);
        }},
        {Generator::Assembler, [&] () {
            return std::make_unique<AssemblerCandidateVariantGenerator>(reference_.get(),
                                                                        kmer_size_.get(),
                                                                        max_variant_size_);
        }},
        {Generator::External, [&] () {
            return std::make_unique<ExternalCandidateVariantGenerator>(variant_source_);
        }},
        {Generator::Online, [&] () {
            return std::make_unique<OnlineCandidateVariantGenerator>(reference_.get(),
                                                                     max_variant_size_);
        }},
        {Generator::Random, [&] () {
            return std::make_unique<RandomCandidateVariantGenerator>(reference_.get());
        }}
    }
    {}
    
    void CandidateGeneratorBuilder::add_generator(const Generator type)
    {
        if (std::find(std::cbegin(generators_), std::cend(generators_), type) == std::cend(generators_)) {
            generators_.push_back(type);
        }
    }
    
    void CandidateGeneratorBuilder::set_reference(const ReferenceGenome& reference)
    {
        reference_ = reference;
    }
    
    void CandidateGeneratorBuilder::set_min_snp_base_quality(const QualityType quality)
    {
        min_snp_base_quality_ = quality;
    }
    
    void CandidateGeneratorBuilder::set_min_supporting_reads(const unsigned num_reads)
    {
        min_supporting_reads_ = num_reads;
    }
    
    void CandidateGeneratorBuilder::set_max_variant_size(const SizeType size)
    {
        max_variant_size_ = size;
    }
    
    void CandidateGeneratorBuilder::set_kmer_size(const unsigned kmer_size)
    {
        kmer_size_ = kmer_size;
    }
    
    void CandidateGeneratorBuilder::set_variant_source(boost::filesystem::path variant_source)
    {
        variant_source_ = std::make_shared<VcfReader>(std::move(variant_source));
    }
    
    void CandidateGeneratorBuilder::set_variant_source(const std::shared_ptr<VcfReader>& variant_source)
    {
        variant_source_ = variant_source;
    }
    
    CandidateVariantGenerator CandidateGeneratorBuilder::build() const
    {
        CandidateVariantGenerator result {};
        
        for (const auto type : generators_) {
            result.register_generator(generator_factory_.at(type)());
        }
        
        return result;
    }
} // namespace Octopus
