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
, assembler_active_region_generator_ {}
, max_read_length_ {}
{}

void ActiveRegionGenerator::add_generator(const std::string& name)
{
    if (is_assembler(name) && !options_.assemble_all) {
        assembler_active_region_generator_ = AssemblerActiveRegionGenerator {reference_};
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

std::vector<GenomicRegion> ActiveRegionGenerator::generate(const GenomicRegion& region, const std::string& generator) const
{
    if (using_assembler() && max_read_length_ > 0) {
        if (!repeats_ || repeats_->request_region != region) {
            repeats_ = RepeatRegions {region, {}, {}, {}};
            const auto tandem_repeats = find_exact_tandem_repeats(reference_, region, max_read_length_ / 2);
            InexactRepeatDefinition minisatellite_def {};
            minisatellite_def.min_exact_repeat_seed_length = 2 * max_read_length_ / 3;
            minisatellite_def.max_seed_join_distance = max_read_length_ / 3;
            minisatellite_def.min_joined_repeat_length = 2 * max_read_length_;
            repeats_->minisatellites = find_repeat_regions(tandem_repeats, region, minisatellite_def);
            InexactRepeatDefinition compound_microsatellite_def {};
            compound_microsatellite_def.max_exact_repeat_seed_period = 6;
            compound_microsatellite_def.min_exact_repeat_seed_length = 4;
            compound_microsatellite_def.max_seed_join_distance = 2;
            compound_microsatellite_def.min_joined_repeat_length = max_read_length_ / 2;
            repeats_->compound_microsatellites = find_repeat_regions(tandem_repeats, region, compound_microsatellite_def);
            if (!repeats_->minisatellites.empty()) {
                std::cout << "minisatellites: ";
                for (const auto& region : repeats_->minisatellites) std::cout << region << " ";
                std::cout << std::endl;
            }
            if (!repeats_->compound_microsatellites.empty()) {
                std::cout << "compound microsatellites: ";
                for (const auto& region : repeats_->compound_microsatellites) std::cout << region << " ";
                std::cout << std::endl;
            }
        }
        if (assembler_active_region_generator_) {
            if (!assembler_active_regions_ || assembler_active_regions_->request_region != region) {
                assembler_active_regions_ = AssemblerActiveRegions {region, {}};
                assembler_active_regions_->active_regions = assembler_active_region_generator_->generate(region);
            }
            if (!repeats_->compound_microsatellites.empty()) {
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
                std::vector<GenomicRegion> cigar_scanner_ignore {};
                cigar_scanner_ignore.reserve(repeats_->minisatellites.size());
                for (const auto& repeat : repeats_->minisatellites) {
                    cigar_scanner_ignore.push_back(expand(repeat, -static_cast<GenomicRegion::Distance>(max_read_length_ / 2)));
                }
                if (!repeats_->assembler_microsatellites.empty()) {
                    cigar_scanner_ignore = merge(std::move(cigar_scanner_ignore), repeats_->assembler_microsatellites);
                }
                return extract_intervening_regions(cigar_scanner_ignore, region);
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
    return assembler_active_region_generator_.is_initialized();
}

} // namespace coretools
} // namespace octopus
