// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_generator.hpp"

#include <algorithm>
#include <cassert>

#include "utils/append.hpp"
#include "utils/parallel_transform.hpp"
namespace octopus { namespace coretools {

VariantGenerator::VariantGenerator()
: debug_log_ {logging::get_debug_log()}
, trace_log_ {logging::get_trace_log()}
, active_region_generator_ {}
{}

VariantGenerator::VariantGenerator(ActiveRegionGenerator region_generator)
: debug_log_ {logging::get_debug_log()}
, trace_log_ {logging::get_trace_log()}
, active_region_generator_ {std::move(region_generator)}
{}

VariantGenerator::VariantGenerator(const VariantGenerator& other)
{
    this->variant_generators_.reserve(other.variant_generators_.size());
    for (const auto& generator : other.variant_generators_) {
        this->variant_generators_.push_back(generator->clone());
    }
    this->active_region_generator_ = other.active_region_generator_;
}

VariantGenerator& VariantGenerator::operator=(VariantGenerator other)
{
    swap(*this, other);
    return *this;
}

VariantGenerator::VariantGenerator(VariantGenerator&& other) : VariantGenerator {}
{
    swap(*this, other);
}

void swap(VariantGenerator& lhs, VariantGenerator& rhs) noexcept
{
    using std::swap;
    swap(lhs.variant_generators_, rhs.variant_generators_);
    swap(lhs.active_region_generator_, rhs.active_region_generator_);
}

void VariantGenerator::add(std::unique_ptr<VariantGenerator> generator)
{
    if (active_region_generator_) active_region_generator_->add_generator(generator->name());
    variant_generators_.push_back(std::move(generator));
}

unsigned VariantGenerator::num_generators() const noexcept
{
    return static_cast<unsigned>(variant_generators_.size());
}

std::unique_ptr<VariantGenerator> VariantGenerator::clone() const
{
    return do_clone();
}

namespace debug {

void log_active_regions(const std::vector<GenomicRegion>& regions, const std::string& generator,
                        boost::optional<logging::DebugLogger>& log)
{
    if (log) {
        auto log_stream = stream(*log);
        log_stream << generator << " active regions: ";
        for (const auto& region : regions) log_stream << region << ' ';
    }
}

void log_candidates(const std::vector<Variant>& candidates, const std::string& generator,
                    boost::optional<logging::DebugLogger>& log)
{
    if (log) {
        auto log_stream = stream(*log);
        if (candidates.empty()) {
            log_stream << "No candidates generated from " << generator << '\n';
        } else {
            log_stream << "Generated " << candidates.size();
            log_stream << " candidate";
            if (candidates.size() > 1) log_stream << "s";
            log_stream << " from " << generator << ":\n";
            for (const auto& c : candidates) log_stream << c << '\n';
        }
    }
}

} // namespace debug

std::vector<Variant> VariantGenerator::generate(const GenomicRegion& region, OptionalThreadPool workers) const
{
    std::vector<Variant> result {};
    const auto generate_helper = [&] (const auto& generator) {
        const auto active_regions = generate_active_regions(region, *generator);
        debug::log_active_regions(active_regions, generator->name(), debug_log_);
        auto result = generator->do_generate(active_regions, workers);
        debug::log_candidates(result, generator->name(), debug_log_);
        assert(std::is_sorted(std::cbegin(result), std::cend(result)));
        return result;
    };
    if (workers && variant_generators_.size() > 1) {
        std::vector<std::vector<Variant>> generator_results(variant_generators_.size());
        using octopus::transform;
        transform(std::cbegin(variant_generators_), std::cend(variant_generators_), 
                  std::begin(generator_results), generate_helper, *workers); 
        for (auto& generator_result : generator_results) {
            utils::append(std::move(generator_result), result);
        }
        std::sort(std::begin(result), std::end(result));
    } else {
        for (const auto& generator : variant_generators_) {
            auto itr = utils::append(generate_helper(generator), result);
            std::inplace_merge(std::begin(result), itr, std::end(result));
        }
    }
    // Each generator is guaranteed to return unique variants, but two generators can still
    // propose the same variants independently.
    remove_duplicates(result);
    return result;
}

bool VariantGenerator::requires_reads() const noexcept
{
    return std::any_of(std::cbegin(variant_generators_), std::cend(variant_generators_),
                       [] (const auto& generator) { return generator->do_requires_reads(); });
}

void VariantGenerator::add_read(const SampleName& sample, const AlignedRead& read)
{
    if (active_region_generator_) active_region_generator_->add_read(sample, read);
    for (auto& generator : variant_generators_) generator->do_add_read(sample, read);
}

void VariantGenerator::add_template(const SampleName& sample, const AlignedTemplate& reads)
{
    if (active_region_generator_) active_region_generator_->add_template(sample, reads);
    for (auto& generator : variant_generators_) generator->do_add_template(sample, reads);
}

void VariantGenerator::do_add_template(const SampleName& sample, const AlignedTemplate& reads)
{
    for (const AlignedRead& read : reads) {
        do_add_read(sample, read);
    }
}

void VariantGenerator::do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last)
{
    std::for_each(first, last, [&] (const auto& read) { do_add_read(sample, read); });
}
void VariantGenerator::do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last)
{
    std::for_each(first, last, [&] (const auto& read) { do_add_read(sample, read); });
}
void VariantGenerator::do_add_reads(const SampleName& sample, TemplateVectorIterator first, TemplateVectorIterator last)
{
    std::for_each(first, last, [&] (const auto& reads) { do_add_template(sample, reads); });
}
void VariantGenerator::do_add_reads(const SampleName& sample, TemplateFlatSetIterator first, TemplateFlatSetIterator last)
{
    std::for_each(first, last, [&] (const auto& reads) { do_add_template(sample, reads); });
}

void VariantGenerator::clear() noexcept
{
    if (active_region_generator_) active_region_generator_->clear();
    for (auto& generator : variant_generators_) generator->do_clear();
}

std::unique_ptr<VariantGenerator> VariantGenerator::do_clone() const
{
    return std::make_unique<VariantGenerator>(*this);
}

VariantGenerator::RegionSet
VariantGenerator::generate_active_regions(const GenomicRegion& region, const VariantGenerator& generator) const
{
    if (active_region_generator_) {
        return active_region_generator_->generate(region, generator.name());
    } else {
        return {region};
    }
}

} // namespace coretools
} // namespace octopus
