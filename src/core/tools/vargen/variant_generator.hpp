// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_generator_hpp
#define variant_generator_hpp

#include <vector>
#include <string>
#include <memory>
#include <iterator>
#include <functional>
#include <cstddef>
#include <type_traits>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "logging/logging.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/variant.hpp"
#include "containers/mappable_flat_multi_set.hpp"

namespace octopus {

class GenomicRegion;

namespace coretools {

class VariantGenerator
{
public:
    VariantGenerator();
    
    VariantGenerator(const VariantGenerator&);
    VariantGenerator& operator=(VariantGenerator);
    VariantGenerator(VariantGenerator&&)            = default;
    VariantGenerator& operator=(VariantGenerator&&) = default;
    
    virtual ~VariantGenerator() = default;
    
    void add(std::unique_ptr<VariantGenerator> generator);
    
    unsigned num_generators() const noexcept;
    
    std::unique_ptr<VariantGenerator> clone() const;
    
    std::vector<Variant> generate(const GenomicRegion& region);
    
    bool requires_reads() const noexcept;
    
    void add_read(const SampleName& sample, const AlignedRead& read);
    
    template <typename InputIt>
    void add_reads(const SampleName& sample, InputIt first, InputIt last);
    
    void clear() noexcept;
    
protected:
    using VectorIterator  = std::vector<AlignedRead>::const_iterator;
    using FlatSetIterator = MappableFlatMultiSet<AlignedRead>::const_iterator;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
    
private:
    std::vector<std::unique_ptr<VariantGenerator>> generators_;
    
    virtual std::unique_ptr<VariantGenerator> do_clone() const;
    
    virtual std::vector<Variant> do_generate_variants(const GenomicRegion& region) { return {}; };
    
    virtual bool do_requires_reads() const noexcept { return false; };
    
    virtual void do_add_read(const SampleName& sample, const AlignedRead& read) {};
    
    // add_reads is not strictly necessary as the effect of calling add_reads must be the same as
    // calling add_read for each read. However, there may be performance benefits
    // to having an add_reads method to avoid many virtual dispatches.
    // Ideally add_reads would be a template to accept any InputIterator, but it is not possible
    // to have template virtual methods. The best solution is therefore to just overload add_reads
    // for common container iterators, more can easily be added if needed.
    virtual void do_add_reads(const SampleName& sample, VectorIterator first, VectorIterator last) {};
    virtual void do_add_reads(const SampleName& sample, FlatSetIterator first, FlatSetIterator last) {};
    
    virtual void do_clear() noexcept {};
    
    virtual std::string name() const { return "VariantGenerator"; }
};

template <typename InputIt>
void VariantGenerator::add_reads(const SampleName& sample, InputIt first, InputIt last)
{
    for (auto& generator : generators_) generator->do_add_reads(sample, first, last);
}

// non-member methods

namespace detail {

template <typename Container, typename G>
void add_reads(const Container& reads, G& generator, std::true_type)
{
    generator.add_reads("octopus-sample", std::cbegin(reads), std::cend(reads));
}

template <typename ReadMap, typename G>
void add_reads(const ReadMap& reads, G& generator, std::false_type)
{
    for (const auto& p : reads) {
        generator.add_reads(p.first, std::cbegin(p.second), std::cend(p.second));
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
