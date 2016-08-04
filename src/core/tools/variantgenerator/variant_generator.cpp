// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_generator.hpp"

#include <algorithm>

#include <basics/genomic_region.hpp>
#include <basics/aligned_read.hpp>

// generators
#include "cigar_scanner.hpp"
#include "local_reassembler.hpp"
#include "vcf_extractor.hpp"
#include "downloader.hpp"
#include "randomiser.hpp"

namespace octopus { namespace coretools {

VariantGenerator::VariantGenerator()
:
debug_log_ {logging::get_debug_log()},
trace_log_ {logging::get_trace_log()}
{}

VariantGenerator::VariantGenerator(const VariantGenerator& other)
{
    generators_.reserve(other.generators_.size());
    for (const auto& generator : other.generators_) {
        generators_.push_back(generator->clone());
    }
}

VariantGenerator& VariantGenerator::operator=(VariantGenerator other)
{
    std::swap(generators_, other.generators_);
    return *this;
}

void VariantGenerator::add(std::unique_ptr<VariantGenerator> generator)
{
    generators_.push_back(std::move(generator));
}

unsigned VariantGenerator::num_generators() const noexcept
{
    return static_cast<unsigned>(generators_.size());
}

std::unique_ptr<VariantGenerator> VariantGenerator::clone() const
{
    return do_clone();
}

namespace debug {
    template <typename S, typename Container>
    void print_generated_candidates(S&& stream, const Container& candidates,
                                    const std::string& generator_name)
    {
        if (candidates.empty()) {
            stream << "No candidates generated from " << generator_name << '\n';
        } else {
            stream << "Generated " << candidates.size();
            stream << " candidate";
            if (candidates.size() > 1) stream << "s";
            stream << " from " << generator_name << ":\n";
            for (const auto& c : candidates) stream << c << '\n';
        }
    }
} // namespace debug

std::vector<Variant> VariantGenerator::generate(const GenomicRegion& region)
{
    std::vector<Variant> result {};
    
    for (auto& generator : generators_) {
        auto generator_result = generator->do_generate_variants(region); // results are sorted
        
        if (debug_log_) {
            debug::print_generated_candidates(stream(*debug_log_), generator_result, generator->name());
        }
        
        auto it = result.insert(std::end(result),
                                std::make_move_iterator(std::begin(generator_result)),
                                std::make_move_iterator(std::end(generator_result)));
        
        std::inplace_merge(std::begin(result), it, std::end(result));
    }
    
    remove_duplicates(result);
    
    return result;
}

bool VariantGenerator::requires_reads() const noexcept
{
    return std::any_of(std::cbegin(generators_), std::cend(generators_),
                       [] (const auto& generator) { return generator->do_requires_reads(); });
}

void VariantGenerator::add_read(const AlignedRead& read)
{
    for (auto& generator : generators_) generator->do_add_read(read);
}

void VariantGenerator::clear() noexcept
{
    for (auto& generator : generators_) generator->do_clear();
}

std::unique_ptr<VariantGenerator> VariantGenerator::do_clone() const
{
    return std::make_unique<VariantGenerator>(*this);
}

// VariantGenerator::Builder

VariantGenerator::Builder& VariantGenerator::Builder::add_generator(const Generator type)
{
    if (std::find(std::cbegin(generators_), std::cend(generators_), type) == std::cend(generators_)) {
        generators_.push_back(type);
    }
    return *this;
}

VariantGenerator::Builder& VariantGenerator::Builder::set_min_base_quality(const BaseQuality min)
{
    parameters_.min_base_quality = min;
    return *this;
}

VariantGenerator::Builder& VariantGenerator::Builder::set_min_supporting_reads(const unsigned num_reads)
{
    parameters_.min_supporting_reads = num_reads;
    return *this;
}

VariantGenerator::Builder& VariantGenerator::Builder::set_max_variant_size(const Variant::RegionType::Size max)
{
    parameters_.max_variant_size = max;
    return *this;
}

VariantGenerator::Builder& VariantGenerator::Builder::add_kmer_size(const unsigned kmer_size)
{
    parameters_.kmer_sizes.push_back(kmer_size);
    return *this;
}

VariantGenerator::Builder& VariantGenerator::Builder::set_assembler_min_base_quality(const BaseQuality min)
{
    parameters_.min_assembler_base_quality = min;
    return *this;
}

VariantGenerator::Builder& VariantGenerator::Builder::set_variant_source(boost::filesystem::path variant_source)
{
    parameters_.variant_source = std::make_shared<VcfReader>(std::move(variant_source));
    return *this;
}

VariantGenerator::Builder& VariantGenerator::Builder::set_variant_source(const std::shared_ptr<const VcfReader>& variant_source)
{
    parameters_.variant_source = variant_source;
    return *this;
}

VariantGenerator VariantGenerator::Builder::build(const ReferenceGenome& reference) const
{
    const auto factory = generate_factory();
    
    VariantGenerator result {};
    
    for (const auto type : generators_) {
        result.add(factory.at(type)(reference));
    }
    
    return result;
}

VariantGenerator::Builder::GeneratorFactoryMap VariantGenerator::Builder::generate_factory() const
{
    return GeneratorFactoryMap {
        {Generator::Alignment, [this] (const ReferenceGenome& reference) {
            CigarScanner::Options options {
                parameters_.min_base_quality,
                parameters_.min_supporting_reads,
                parameters_.max_variant_size
            };
            return std::make_unique<CigarScanner>(reference, options);
        }},
        {Generator::Assembler, [this] (const ReferenceGenome& reference) {
            auto quality = (parameters_.min_assembler_base_quality)
            ? *parameters_.min_assembler_base_quality : parameters_.min_base_quality;
            LocalReassembler::Options options {
                parameters_.kmer_sizes,
                quality,
                parameters_.min_supporting_reads,
                parameters_.max_variant_size
            };
            return std::make_unique<LocalReassembler>(reference, std::move(options));
        }},
        {Generator::External, [this] (const ReferenceGenome& reference) {
            return std::make_unique<VcfExtractor>(parameters_.variant_source);
        }},
        {Generator::Online, [this] (const ReferenceGenome& reference) {
            return std::make_unique<Downloader>(reference, parameters_.max_variant_size);
        }},
        {Generator::Random, [this] (const ReferenceGenome& reference) {
            return std::make_unique<Randomiser>(reference);
        }}
    };
}

} // namespace coretools
} // namespace octopus
