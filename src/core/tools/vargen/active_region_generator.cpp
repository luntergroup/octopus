// Copyright (c) 2017 Daniel Cooke
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
            assembler_active_region_generator_ = AssemblerActiveRegionGenerator {reference_};
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
    if (using_assembler() && max_read_length_ > 0) {
        if (!repeats_ || repeats_->request_region != region) {
            repeats_ = RepeatRegions {region, {}, {}, {}};
            const auto tandem_repeats = find_exact_tandem_repeats(reference_, region, max_read_length_ / 2);
            repeats_->minisatellites = find_minisatellites(tandem_repeats, region, max_read_length_);
            repeats_->compound_microsatellites = find_compound_microsatellites(tandem_repeats, region, max_read_length_);
        }
        if (assembler_active_region_generator_) {
            if (!assembler_active_regions_ || assembler_active_regions_->request_region != region) {
                assembler_active_regions_ = AssemblerActiveRegions {region, {}};
                assembler_active_regions_->active_regions = assembler_active_region_generator_->generate(region);
            }
            if (!repeats_->compound_microsatellites.empty() && repeats_->assembler_microsatellites.empty()) {
                repeats_->assembler_microsatellites.reserve(repeats_->compound_microsatellites.size());
                auto microsatellites_begin_itr = std::cbegin(repeats_->compound_microsatellites);
                const auto microsatellites_end_itr = std::cend(repeats_->compound_microsatellites);
                for (const auto& assembler_region : assembler_active_regions_->active_regions) {
                    if (microsatellites_begin_itr == microsatellites_end_itr) break;
                    const auto overlapped_microsatellites = bases(overlap_range(microsatellites_begin_itr, microsatellites_end_itr,
                                                                                assembler_region, BidirectionallySortedTag {}));
                    append(overlapped_microsatellites, repeats_->assembler_microsatellites);
                    microsatellites_begin_itr = overlapped_microsatellites.end();
                }
            }
        }
        if (is_assembler(generator)) {
            if (assembler_active_region_generator_) {
                auto result = merge(std::move(assembler_active_regions_->active_regions), repeats_->minisatellites);
                assembler_active_regions_ = boost::none;
                if (repeats_->assembler_microsatellites.empty()) {
                    return result;
                } else {
                    return merge(std::move(result), repeats_->assembler_microsatellites);
                }
            } else {
                return {region};
            }
        } else if (is_cigar_scanner(generator)) {
            if (repeats_->minisatellites.empty() && repeats_->assembler_microsatellites.empty()) {
                return {region};
            } else {
                std::vector<GenomicRegion> repeat_regions {};
                repeat_regions.reserve(repeats_->minisatellites.size() + repeats_->assembler_microsatellites.size());
                const auto max_read_distance = static_cast<GenomicRegion::Distance>(max_read_length_);
                for (const auto& repeat : repeats_->minisatellites) {
                    if (size(repeat) > 3 * max_read_length_) {
                        repeat_regions.push_back(expand(repeat, -max_read_distance));
                    } else {
                        assert(size(repeat) > 2 * max_read_length_);
                        repeat_regions.push_back(expand(repeat, -max_read_distance / 2));
                    }
                }
                for (const auto& repeat : repeats_->assembler_microsatellites) {
                    if (size(repeat) > 3 * max_read_length_) {
                        repeat_regions.push_back(expand(repeat, -max_read_distance));
                    } else if (size(repeat) > 2 * max_read_length_) {
                        repeat_regions.push_back(expand(repeat, -max_read_distance / 2));
                    } else {
                        repeat_regions.push_back(repeat);
                    }
                }
                if (!repeat_regions.empty()) {
                    std::sort(std::begin(repeat_regions), std::end(repeat_regions));
                    return extract_intervening_regions(extract_covered_regions(repeat_regions), region);
                } else {
                    return {region};
                }
            }
        }
    }
    return {region};
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
