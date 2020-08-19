// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "active_region_generator.hpp"

#include <iterator>

#include "utils/repeat_finder.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/append.hpp"

namespace octopus { namespace coretools {

ActiveRegionGenerator::ActiveRegionGenerator(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
, assembler_name_ {"LocalReassembler"}
, cigar_scanner_name_ {"CigarScanner"}
, using_assembler_ {false}
, assembler_active_region_generator_ {}
, max_read_length_ {}
{}

void ActiveRegionGenerator::add_generator(const std::string& name)
{
    if (is_assembler(name)) {
        using_assembler_ = true;
        if (!options_.assemble_all) {
            if (options_.assembler_active_region_generator_options) {
                assembler_active_region_generator_ = AssemblerActiveRegionGenerator {reference_, *options_.assembler_active_region_generator_options};
            } else {
                assembler_active_region_generator_ = AssemblerActiveRegionGenerator {reference_};
            }
        }
    }
}

void ActiveRegionGenerator::add_read(const SampleName& sample, const AlignedRead& read)
{
    max_read_length_ = std::max(max_read_length_, sequence_size(read));
    if (assembler_active_region_generator_) assembler_active_region_generator_->add(sample, read);
}

auto merge(std::vector<GenomicRegion> lhs, std::vector<GenomicRegion> rhs)
{
    auto itr = utils::append(std::move(rhs), lhs);
    std::inplace_merge(std::begin(lhs), itr, std::end(lhs));
    return extract_covered_regions(lhs);
}

template <typename Range, typename T>
auto append(const Range& range, std::vector<T>& result)
{
    return result.insert(std::cend(result), std::cbegin(range), std::cend(range));
};

auto find_minisatellites(const std::vector<TandemRepeat>& repeats, const GenomicRegion& region,
                         const std::size_t max_read_length)
{
    InexactRepeatDefinition minisatellite_def {};
    minisatellite_def.min_exact_repeat_seed_length = 2 * max_read_length / 3;
    minisatellite_def.min_exact_repeat_seed_periods = 3;
    minisatellite_def.max_seed_join_distance = max_read_length / 3;
    minisatellite_def.min_joined_repeat_length = 2 * max_read_length;
    return join(find_repeat_regions(repeats, region, minisatellite_def), max_read_length / 2);
}

auto find_compound_microsatellites(const std::vector<TandemRepeat>& repeats, const GenomicRegion& region,
                                   const std::size_t max_read_length)
{
    InexactRepeatDefinition compound_microsatellite_def {};
    compound_microsatellite_def.max_exact_repeat_seed_period = 6;
    compound_microsatellite_def.min_exact_repeat_seed_length = 4;
    compound_microsatellite_def.min_exact_repeat_seed_periods = 4;
    compound_microsatellite_def.min_exact_seeds = 2;
    compound_microsatellite_def.max_seed_join_distance = 1;
    compound_microsatellite_def.min_joined_repeat_length = max_read_length / 4;
    return join(find_repeat_regions(repeats, region, compound_microsatellite_def), max_read_length / 2);
}

std::vector<GenomicRegion> ActiveRegionGenerator::generate(const GenomicRegion& region, const std::string& generator) const
{
    if (is_assembler(generator) && assembler_active_region_generator_) {
        return assembler_active_region_generator_->generate(region);
    } else {
        return {region};
    }
}

void ActiveRegionGenerator::clear() noexcept
{
    if (assembler_active_region_generator_) assembler_active_region_generator_->clear();
}

// private methods

bool ActiveRegionGenerator::is_cigar_scanner(const std::string& generator) const noexcept
{
    return generator == cigar_scanner_name_;
}

bool ActiveRegionGenerator::is_assembler(const std::string& generator) const noexcept
{
    return generator == assembler_name_;
}

bool ActiveRegionGenerator::using_assembler() const noexcept
{
    return using_assembler_;
}

} // namespace coretools
} // namespace octopus
