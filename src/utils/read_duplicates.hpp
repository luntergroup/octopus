// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_duplicates_hpp
#define read_duplicates_hpp

#include <iterator>
#include <algorithm>
#include <vector>
#include <map>

#include "basics/aligned_read.hpp"
#include "basics/aligned_template.hpp"
#include "utils/append.hpp"

namespace octopus {

struct FivePrimeDuplicateDefinition
{
    bool unpaired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
    bool paired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
};

struct FivePrimeAndCigarDuplicateDefinition
{
    bool unpaired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
    bool paired_equal(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
};

namespace detail {
	
template <typename Range, typename DuplicateDefinition>
auto find_duplicate(const AlignedRead& read, const Range& reads, const DuplicateDefinition& dup_def) noexcept
{
	return std::find_if(std::begin(reads), std::end(reads), [&] (auto other_itr) { return dup_def.paired_equal(read, *other_itr); });
}

struct NextSegmentLess
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
};

} // namespace detail

template <typename ForwardIt,
          typename DuplicateDefinition>
std::vector<std::vector<ForwardIt>>
find_duplicate_reads(ForwardIt first, const ForwardIt last,
                     const DuplicateDefinition& duplicate_definition)
{
    std::vector<std::vector<ForwardIt>> result {};
    std::vector<ForwardIt> buffer {};
    // Recall that reads come sorted w.r.t operator< and it is therefore not guaranteed that 'duplicate'
    // reads (according to IsDuplicate) will be adjacent to one another. In particular, operator< only
    // guarantees that duplicate read segment described in the AlignedRead object will be adjacent.
    const auto are_primary_dups = [&] (const auto& lhs, const auto& rhs) { return duplicate_definition.unpaired_equal(lhs, rhs); };
    for (first = std::adjacent_find(first, last, are_primary_dups); first != last; first = std::adjacent_find(first, last, are_primary_dups)) {
        buffer.reserve(8);
        if (first->has_other_segment()) buffer.push_back(first);
        if (std::next(first)->has_other_segment()) buffer.push_back(std::next(first));
        const AlignedRead& primary_dup_read {*first};
        std::advance(first, 2);
        for (; first != last && are_primary_dups(*first, primary_dup_read); ++first) {
            if (first->has_other_segment()) buffer.push_back(first);
        }
        std::sort(buffer.begin(), buffer.end(), [] (ForwardIt lhs, ForwardIt rhs) { return detail::NextSegmentLess{}(*lhs, *rhs); });
        const auto are_duplicates = [&] (ForwardIt lhs, ForwardIt rhs) { return duplicate_definition.paired_equal(*lhs, *rhs); };
        for (auto dup_itr = std::adjacent_find(buffer.cbegin(), buffer.cend(), are_duplicates); dup_itr != buffer.cend();
                 dup_itr = std::adjacent_find(dup_itr, buffer.cend(), are_duplicates)) {
            std::vector<ForwardIt> duplicates {dup_itr, std::next(dup_itr, 2)};
            const AlignedRead& duplicate_read {**dup_itr};
            std::advance(dup_itr, 2);
            for (; dup_itr != buffer.cend() && duplicate_definition.paired_equal(**dup_itr, duplicate_read); ++dup_itr) {
                duplicates.push_back(*dup_itr);
            }
            result.push_back(std::move(duplicates));
        }
        buffer.clear();
    }
    return result;
}

template <typename ForwardIt>
std::vector<std::vector<ForwardIt>>
find_duplicate_reads(ForwardIt first, const ForwardIt last)
{
    return find_duplicate_reads(first, last, FivePrimeDuplicateDefinition {});
}

namespace detail {

template <typename AlignedReadIterator>
struct AlignedReadIteratorNameLess
{
	using is_transparent = void;
	bool operator()(const AlignedReadIterator& lhs, const AlignedReadIterator& rhs) const noexcept { return lhs->name() < rhs->name(); }
	bool operator()(const AlignedReadIterator& lhs, const AlignedRead& rhs) const noexcept { return lhs->name() < rhs.name(); }
	bool operator()(const AlignedRead& lhs, const AlignedReadIterator& rhs) const noexcept { return lhs.name() < rhs->name(); }
	bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept { return lhs.name() < rhs.name(); }
};

} // namespace detail

template <typename ForwardIt,
          typename BinaryPredicate,
          typename DuplicateDefinition>
ForwardIt 
remove_duplicate_reads(ForwardIt first, ForwardIt last,
                       const BinaryPredicate duplicate_compare,
                       const DuplicateDefinition& duplicate_definition)
{
    // See comment in 'find_duplicates'
    const auto are_primary_dups = [&] (const auto& lhs, const auto& rhs) { return duplicate_definition.unpaired_equal(lhs, rhs); };
    first = std::adjacent_find(first, last, are_primary_dups);
    if (first != last) {
        std::vector<ForwardIt> candidate_duplicates {first++};
		using DuplicatePairedReadMap = std::map<ForwardIt, std::vector<AlignedRead>, detail::AlignedReadIteratorNameLess<ForwardIt>>;
		DuplicatePairedReadMap paired_duplicates {}, working_paired_duplicates {};
        candidate_duplicates.reserve(100);
        for (auto read_itr = first; read_itr != last; ++read_itr) {
			AlignedRead& read {*read_itr};
            if (are_primary_dups(read, *candidate_duplicates.front())) { // can check any candidate
				const auto duplicate_itr = detail::find_duplicate(read, candidate_duplicates, duplicate_definition);
                if (duplicate_itr == std::end(candidate_duplicates)) { // read may not be a duplicate
                    if (read_itr != first) *first = std::move(read);
                    candidate_duplicates.emplace_back(first++);
                } else { // read is a duplicate
					const ForwardIt curr_best_duplicate_itr {*duplicate_itr};
					if (duplicate_compare(*curr_best_duplicate_itr, read)) {
						// swap duplicates
						if (read.has_other_segment()) {
							assert(curr_best_duplicate_itr->has_other_segment());
							const auto replace_itr = working_paired_duplicates.find(curr_best_duplicate_itr);
							if (replace_itr != std::cend(working_paired_duplicates)) {
								assert(replace_itr->first == curr_best_duplicate_itr);
								auto worse_duplicates = std::move(replace_itr->second);
								worse_duplicates.push_back(std::move(*replace_itr->first));
								working_paired_duplicates.erase(replace_itr);
								*curr_best_duplicate_itr = std::move(read);
								working_paired_duplicates[curr_best_duplicate_itr] = std::move(worse_duplicates);
							} else {
								std::swap(*curr_best_duplicate_itr, read);
								working_paired_duplicates[curr_best_duplicate_itr] = {std::move(read)};
							}
						} else {
							std::iter_swap(*duplicate_itr, read_itr);
						}
					} else if (read.has_other_segment()) {
						working_paired_duplicates[curr_best_duplicate_itr].push_back(std::move(read));
					}
                }
            } else {
                if (read_itr != first) *first = std::move(read);
                candidate_duplicates.assign({first++});
				for (auto& p : working_paired_duplicates) {
					auto mate_itr = paired_duplicates.find(*p.first);
					if (mate_itr != std::cend(paired_duplicates)) {
						// read's mate was the best duplicate segment and so is read - all is good
						paired_duplicates.erase(mate_itr);
					} else {
						// Either read's mate (and its duplicates) were not considered (maybe it was already filrered),
						// or one of read's mate duplicates was 'better'.
						// To find out we need to check all of read's duplicate segments.
						for (auto& duplicate : p.second) {
							mate_itr = paired_duplicates.find(duplicate);
							if (mate_itr != std::cend(paired_duplicates)) break;
						}
						if (mate_itr == std::cend(paired_duplicates)) {
							paired_duplicates.emplace(p.first, std::move(p.second));
						} else {
							// Unpaired duplicate segments were chosen, we need to pick the best pair
							std::vector<AlignedRead> buffer {};
							buffer.reserve(mate_itr->second.size() + p.second.size() + 2);
							buffer.push_back(std::move(*mate_itr->first));
							utils::append(std::move(mate_itr->second), buffer);
							buffer.push_back(std::move(*p.first));
							utils::append(std::move(p.second), buffer);
							std::vector<AlignedTemplate> duplicate_pairs {};
							duplicate_pairs.reserve(buffer.size() / 2);
							ReadLinkageConfig template_config {};
        					template_config.linkage = ReadLinkageType::paired;
							template_config.linked_only = true;
							make_read_templates(std::cbegin(buffer), std::cend(buffer), std::back_inserter(duplicate_pairs), template_config);
							assert(duplicate_pairs.size() > 0);
							const auto& best_pair = *std::max_element(std::cbegin(duplicate_pairs), std::cend(duplicate_pairs), duplicate_compare);
							assert(best_pair.size() == 2);
							*mate_itr->first = best_pair[0];
							*p.first = best_pair[1];
							paired_duplicates.erase(mate_itr);
						}
					}
				}
				working_paired_duplicates.clear();
            }
        }
    }
    return first;
}

struct DuplicateReadLess
{
	bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
	{
		return lhs.mapping_quality() < rhs.mapping_quality()
			 || (lhs.mapping_quality() == rhs.mapping_quality() && sum_base_qualities(lhs) < sum_base_qualities(rhs)); 
	}
	bool operator()(const AlignedTemplate& lhs, const AlignedTemplate& rhs) const noexcept
	{
		const auto lhs_total_mq = sum_mapping_qualities(lhs);
		const auto rhs_total_mq = sum_mapping_qualities(rhs);
		return lhs_total_mq < rhs_total_mq
			 || (lhs_total_mq == rhs_total_mq && sum_base_qualities(lhs) < sum_base_qualities(rhs)); 
	}
};

template <typename ForwardIt,
          typename DuplicateDefinition>
ForwardIt 
remove_duplicate_reads(ForwardIt first, ForwardIt last,
                       const DuplicateDefinition& duplicate_definition)
{
	return remove_duplicate_reads(first, last, DuplicateReadLess {}, duplicate_definition);
}

template <typename ForwardIt>
ForwardIt
remove_duplicate_reads(ForwardIt first, ForwardIt last)
{
    return remove_duplicate_reads(first, last, FivePrimeDuplicateDefinition {});
}

} // namespace

#endif //read_duplicates_hpp
