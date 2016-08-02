//
//  variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_generator__
#define Octopus_variant_generator__

#include <vector>
#include <string>
#include <deque>
#include <unordered_map>
#include <memory>
#include <iterator>
#include <functional>
#include <cstddef>
#include <type_traits>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "common.hpp"
#include "logging.hpp"
#include "variant.hpp"
#include "mappable_flat_multi_set.hpp"
#include "reference_genome.hpp"
#include "vcf_reader.hpp"

namespace octopus {

class AlignedRead;
class GenomicRegion;

namespace coretools {

class VariantGenerator
{
public:
    class Builder;
    
    VariantGenerator();
    
    VariantGenerator(const VariantGenerator&);
    VariantGenerator& operator=(VariantGenerator);
    VariantGenerator(VariantGenerator&&)                 = default;
    VariantGenerator& operator=(VariantGenerator&&)      = default;
    
    virtual ~VariantGenerator() = default;
    
    void add(std::unique_ptr<VariantGenerator> generator);
    
    unsigned num_generators() const noexcept;
    
    std::unique_ptr<VariantGenerator> clone() const;
    
    std::vector<Variant> generate(const GenomicRegion& region);
    
    bool requires_reads() const noexcept;
    
    void add_read(const AlignedRead& read);
    
    template <typename InputIt>
    void add_reads(InputIt first, InputIt last);
    
    void clear() noexcept;
    
protected:
    using VectorIterator  = std::vector<AlignedRead>::const_iterator;
    using FlatSetIterator = MappableFlatMultiSet<AlignedRead>::const_iterator;
    
    boost::optional<logging::DebugLogger> debug_log_;
    boost::optional<logging::TraceLogger> trace_log_;
    
private:
    std::vector<std::unique_ptr<VariantGenerator>> generators_;
    
    virtual std::unique_ptr<VariantGenerator> do_clone() const;
    
    virtual std::vector<Variant> do_generate_variants(const GenomicRegion& region) { return {}; };
    
    virtual bool do_requires_reads() const noexcept { return false; };
    
    virtual void do_add_read(const AlignedRead& read) {};
    
    // add_reads is not strictly necessary as the effect of calling add_reads must be the same as
    // calling add_read for each read. However, there may be performance benifits
    // to having an add_reads method to avoid many virtual dispatches.
    // Ideally add_reads would be templated to accept any InputIterator, but it is not possible
    // to have template virtual methods. The best solution is therefore to just overload add_reads
    // for common container iterators, more can easily be added if needed.
    virtual void do_add_reads(VectorIterator first, VectorIterator last) {};
    virtual void do_add_reads(FlatSetIterator first, FlatSetIterator last) {};
    
    virtual void do_clear() noexcept {};
    
    virtual std::string name() const { return "Composer"; }
};

class VariantGenerator::Builder
{
public:
    using BaseQuality = AlignedRead::BaseQuality;
    using Position    = GenomicRegion::Position;
    
    enum class Generator { Alignment, Assembler, External, Online, Random };
    
    Builder() = default;
    
    Builder(const Builder&)            = default;
    Builder& operator=(const Builder&) = default;
    Builder(Builder&&)                 = default;
    Builder& operator=(Builder&&)      = default;
    
    ~Builder() = default;
    
    Builder& add_generator(Generator type);
    
    Builder& set_min_base_quality(BaseQuality min);
    Builder& set_min_supporting_reads(unsigned num_reads);
    Builder& set_max_variant_size(Position size);
    
    Builder& add_kmer_size(unsigned kmer_size);
    Builder& set_assembler_min_base_quality(BaseQuality min);
    
    Builder& set_variant_source(boost::filesystem::path source);
    Builder& set_variant_source(const std::shared_ptr<const VcfReader>& source);
    
    VariantGenerator build(const ReferenceGenome& reference) const;
    
private:
    struct Parameters
    {
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
    
    using GeneratorFactory    = std::function<std::unique_ptr<VariantGenerator>(const ReferenceGenome&)>;
    using GeneratorFactoryMap = std::unordered_map<Generator, GeneratorFactory, GeneratorTypeHash>;
    
    std::deque<Generator> generators_;
    
    Parameters parameters_;
    
    GeneratorFactoryMap generate_factory() const;
};

template <typename InputIt>
void VariantGenerator::add_reads(InputIt first, InputIt last)
{
    for (auto& generator : generators_) generator->add_reads(first, last);
}

// non-member methods

namespace detail {
    template <typename Container, typename G>
    void add_reads(const Container& reads, G& generator, std::true_type)
    {
        generator.add_reads(std::cbegin(reads), std::cend(reads));
    }
    
    template <typename ReadMap, typename G>
    void add_reads(const ReadMap& reads, G& generator, std::false_type)
    {
        for (const auto& p : reads) {
            generator.add_reads(std::cbegin(p.second), std::cend(p.second));
        }
    }
} // namespace detail

template <typename Container, typename G,
          typename = std::enable_if_t<std::is_base_of<VariantGenerator, G>::value>>
void add_reads(const Container& reads, G& generator)
{
    detail::add_reads(reads, generator, std::is_same<typename Container::value_type, AlignedRead> {});
}

} // namespace coretools

using coretools::VariantGenerator;

} // namespace octopus

#endif
