// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_generator.hpp"

#include <algorithm>
#include <cassert>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"

namespace octopus { namespace coretools {

VariantGenerator::VariantGenerator()
: debug_log_ {logging::get_debug_log()}
, trace_log_ {logging::get_trace_log()}
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
        auto generator_result = generator->do_generate_variants(region);
        assert(std::is_sorted(std::cbegin(generator_result), std::cend(generator_result)));
        if (debug_log_) {
            debug::print_generated_candidates(stream(*debug_log_), generator_result, generator->name());
        }
        auto itr = result.insert(std::end(result),
                                 std::make_move_iterator(std::begin(generator_result)),
                                 std::make_move_iterator(std::end(generator_result)));
        std::inplace_merge(std::begin(result), itr, std::end(result));
    }
    // Each generator is guaranteed to return unique variants, but two generators can still
    // propose the same variants independently.
    remove_duplicates(result);
    return result;
}

bool VariantGenerator::requires_reads() const noexcept
{
    return std::any_of(std::cbegin(generators_), std::cend(generators_),
                       [] (const auto& generator) { return generator->do_requires_reads(); });
}

void VariantGenerator::add_read(const SampleName& sample, const AlignedRead& read)
{
    for (auto& generator : generators_) generator->do_add_read(sample, read);
}

void VariantGenerator::clear() noexcept
{
    for (auto& generator : generators_) generator->do_clear();
}

std::unique_ptr<VariantGenerator> VariantGenerator::do_clone() const
{
    return std::make_unique<VariantGenerator>(*this);
}


} // namespace coretools
} // namespace octopus
