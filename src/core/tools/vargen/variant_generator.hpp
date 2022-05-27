// Copyright (c) 2015-2021 Daniel Cooke
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
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "basics/aligned_template.hpp"
#include "core/types/variant.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "utils/type_tricks.hpp"
#include "utils/thread_pool.hpp"
#include "active_region_generator.hpp"

namespace octopus { namespace coretools {
class VariantGenerator
{
public:
    using OptionalThreadPool = boost::optional<ThreadPool&>;

    VariantGenerator();
    VariantGenerator(ActiveRegionGenerator region_generator);
    
    VariantGenerator(const VariantGenerator&);
    VariantGenerator& operator=(VariantGenerator);
    VariantGenerator(VariantGenerator&&);
    
    virtual ~VariantGenerator() = default;
    
    friend void swap(VariantGenerator& lhs, VariantGenerator& rhs) noexcept;
    
    void add(std::unique_ptr<VariantGenerator> generator);
    
    unsigned num_generators() const noexcept;
    
    std::unique_ptr<VariantGenerator> clone() const;
    
    std::vector<Variant> generate(const GenomicRegion& region, OptionalThreadPool = boost::none) const;
    
    bool requires_reads() const noexcept;
    
    void add_read(const SampleName& sample, const AlignedRead& read);
    void add_template(const SampleName& sample, const AlignedTemplate& reads);
    template <typename InputIt>
    void add_reads(const SampleName& sample, InputIt first, InputIt last);
    
    void clear() noexcept;
    
protected:
    using ReadVectorIterator  = std::vector<AlignedRead>::const_iterator;
    using ReadFlatSetIterator = MappableFlatMultiSet<AlignedRead>::const_iterator;
    using TemplateVectorIterator = std::vector<AlignedTemplate>::const_iterator;
    using TemplateFlatSetIterator = MappableFlatMultiSet<AlignedTemplate>::const_iterator;
    
    using RegionSet = std::vector<GenomicRegion>;
    
    mutable boost::optional<logging::DebugLogger> debug_log_;
    mutable boost::optional<logging::TraceLogger> trace_log_;
    
private:
    std::vector<std::unique_ptr<VariantGenerator>> variant_generators_;
    boost::optional<ActiveRegionGenerator> active_region_generator_;
    
    virtual std::unique_ptr<VariantGenerator> do_clone() const;
    virtual std::vector<Variant> do_generate(const RegionSet& regions, OptionalThreadPool workers) const { return {}; };
    virtual bool do_requires_reads() const noexcept { return false; };
    virtual void do_add_read(const SampleName& sample, const AlignedRead& read) {};
    virtual void do_add_template(const SampleName& sample, const AlignedTemplate& reads);
    // add_reads is not strictly necessary as the effect of calling add_reads must be the same as
    // calling add_read for each read. However, there may be performance benefits
    // to having an add_reads method to avoid many virtual dispatches.
    // Ideally add_reads would be a template to accept any InputIterator, but it is not possible
    // to have template virtual methods. The best solution is therefore to just overload add_reads
    // for common container iterators, more can easily be added if needed.
    virtual void do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last);
    virtual void do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last);
    virtual void do_add_reads(const SampleName& sample, TemplateVectorIterator first, TemplateVectorIterator last);
    virtual void do_add_reads(const SampleName& sample, TemplateFlatSetIterator first, TemplateFlatSetIterator last);
    virtual void do_clear() noexcept {};
    
    virtual std::string name() const { return "VariantGenerator"; }
    
    RegionSet generate_active_regions(const GenomicRegion& region, const VariantGenerator& generator) const;
};

template <typename InputIt>
void VariantGenerator::add_reads(const SampleName& sample, InputIt first, InputIt last)
{
    if (active_region_generator_) active_region_generator_->add_reads(sample, first, last);
    for (auto& generator : variant_generators_) generator->do_add_reads(sample, first, last);
}

// non-member methods

namespace detail {

template <typename Container, typename G>
void add_reads(const Container& reads, G& generator, std::true_type)
{
    for (const auto& p : reads) {
        generator.add_reads(p.first, std::cbegin(p.second), std::cend(p.second));
    }
}

template <typename ReadMap, typename G>
void add_reads(const ReadMap& reads, G& generator, std::false_type)
{
    generator.add_reads("octopus", std::cbegin(reads), std::cend(reads));
}

} // namespace detail

template <typename Container, typename G,
          typename = std::enable_if_t<std::is_base_of<VariantGenerator, G>::value>>
void add_reads(const Container& reads, G& generator)
{
    detail::add_reads(reads, generator, IsMap<Container> {});
}

} // namespace coretools

using coretools::VariantGenerator;

} // namespace octopus

#endif
