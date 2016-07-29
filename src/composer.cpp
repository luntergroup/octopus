//
//  composer.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "composer.hpp"

#include <algorithm>
#include <iterator>

#include "genomic_region.hpp"
#include "aligned_read.hpp"

#include "cigar_scanner.hpp"
#include "local_reassembler.hpp"
#include "source_file_reader.hpp"
#include "online_crawler.hpp"
#include "randomiser.hpp"

namespace octopus { namespace core { namespace generators
{
void Composer::register_generator(std::unique_ptr<CandidateVariantGenerator> generator)
{
    generators_.emplace_back(std::move(generator));
}

bool Composer::requires_reads() const noexcept
{
    return std::any_of(std::cbegin(generators_), std::cend(generators_),
                       [] (const auto& generator) { return generator->requires_reads(); });
}

void Composer::add_read(const AlignedRead& read)
{
    for (auto& generator : generators_) {
        generator->add_read(read);
    }
}

void Composer::add_reads(std::vector<AlignedRead>::const_iterator first,
                                          std::vector<AlignedRead>::const_iterator last)
{
    for (auto& generator : generators_) {
        generator->add_reads(first, last);
    }
}

void Composer::add_reads(MappableFlatMultiSet<AlignedRead>::const_iterator first,
                                          MappableFlatMultiSet<AlignedRead>::const_iterator last)
{
    for (auto& generator : generators_) {
        generator->add_reads(first, last);
    }
}

std::vector<Variant> Composer::generate_candidates(const GenomicRegion& region)
{
    std::vector<Variant> result {};
    
    for (auto& generator : generators_) {
        auto generator_result = generator->generate_candidates(region); // results are sorted
        auto it = result.insert(std::end(result),
                                std::make_move_iterator(std::begin(generator_result)),
                                std::make_move_iterator(std::end(generator_result)));
        std::inplace_merge(std::begin(result), it, std::end(result));
    }
    
    remove_duplicates(result);
    
    return result;
}

void Composer::reserve(const std::size_t n)
{
    for (auto& generator : generators_) {
        generator->reserve(n);
    }
}

void Composer::clear()
{
    for (auto& generator : generators_) {
        generator->clear();
    }
}

// Composer::Builder

Composer::Builder::Builder()
:
generators_ {},
parameters_ {},
generator_factory_ {generate_factory()}
{}

Composer::Builder::Builder(const Composer::Builder& other)
:
generators_        {other.generators_},
parameters_        {other.parameters_},
generator_factory_ {generate_factory()}
{}

Composer::Builder& Composer::Builder::operator=(const Composer::Builder& other)
{
    generators_        = other.generators_;
    parameters_        = other.parameters_;
    generator_factory_ = generate_factory();
    return *this;
}

Composer::Builder::Builder(Composer::Builder&& other)
:
generators_        {std::move(other.generators_)},
parameters_        {std::move(other.parameters_)},
generator_factory_ {generate_factory()}
{}

Composer::Builder& Composer::Builder::operator=(Composer::Builder&& other)
{
    using std::swap;
    swap(generators_, other.generators_);
    swap(parameters_, other.parameters_);
    generator_factory_ = generate_factory();
    return *this;
}

unsigned Composer::Builder::num_generators() const noexcept
{
    return static_cast<unsigned>(generators_.size());
}

Composer::Builder& Composer::Builder::add_generator(const Generator type)
{
    if (std::find(std::cbegin(generators_), std::cend(generators_), type) == std::cend(generators_)) {
        generators_.push_back(type);
    }
    return *this;
}

Composer::Builder& Composer::Builder::set_reference(const ReferenceGenome& reference)
{
    parameters_.reference = reference;
    return *this;
}

Composer::Builder& Composer::Builder::set_min_base_quality(const BaseQuality min)
{
    parameters_.min_base_quality = min;
    return *this;
}

Composer::Builder& Composer::Builder::set_min_supporting_reads(const unsigned num_reads)
{
    parameters_.min_supporting_reads = num_reads;
    return *this;
}

Composer::Builder& Composer::Builder::set_max_variant_size(const Variant::RegionType::Size max)
{
    parameters_.max_variant_size = max;
    return *this;
}

Composer::Builder& Composer::Builder::add_kmer_size(const unsigned kmer_size)
{
    parameters_.kmer_sizes.push_back(kmer_size);
    return *this;
}

Composer::Builder& Composer::Builder::set_assembler_min_base_quality(const BaseQuality min)
{
    parameters_.min_assembler_base_quality = min;
    return *this;
}

Composer::Builder& Composer::Builder::set_variant_source(boost::filesystem::path variant_source)
{
    parameters_.variant_source = std::make_shared<VcfReader>(std::move(variant_source));
    return *this;
}

Composer::Builder& Composer::Builder::set_variant_source(const std::shared_ptr<const VcfReader>& variant_source)
{
    parameters_.variant_source = variant_source;
    return *this;
}

Composer Composer::Builder::build() const
{
    Composer result {};
    
    for (const auto type : generators_) {
        result.register_generator(generator_factory_.at(type)());
    }
    
    return result;
}

// private methods

Composer::Builder::GeneratorFactoryMap Composer::Builder::generate_factory() const
{
    return GeneratorFactoryMap {
        {Generator::Alignment, [this] () {
            CigarScanner::Options options {
                parameters_.min_base_quality,
                parameters_.min_supporting_reads,
                parameters_.max_variant_size
            };
            return std::make_unique<CigarScanner>(*parameters_.reference, options);
        }},
        {Generator::Assembler, [this] () {
            const auto quality = (parameters_.min_assembler_base_quality)
            ? *parameters_.min_assembler_base_quality : parameters_.min_base_quality;
            LocalReassembler::Options options {
                parameters_.kmer_sizes,
                quality,
                parameters_.min_supporting_reads,
                parameters_.max_variant_size
            };
            return std::make_unique<LocalReassembler>(*parameters_.reference, std::move(options));
        }},
        {Generator::External, [this] () {
            return std::make_unique<SourceFileReader>(parameters_.variant_source);
        }},
        {Generator::Online, [this] () {
            return std::make_unique<OnlineCrawler>(*parameters_.reference,
                                                   parameters_.max_variant_size);
        }},
        {Generator::Random, [this] () {
            return std::make_unique<Randomiser>(*parameters_.reference);
        }}
    };
}
} // namespace generators
} // namespace core
} // namespace octopus
