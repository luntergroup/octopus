// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_duplicates_hpp
#define read_duplicates_hpp

#include <iterator>
#include <algorithm>
#include <vector>

#include "basics/aligned_read.hpp"

namespace octopus {

bool primary_segments_are_duplicates(const AlignedRead& lhs, const AlignedRead& rhs) noexcept;
bool other_segments_are_duplicates(const AlignedRead& lhs, const AlignedRead& rhs) noexcept;
bool are_duplicates(const AlignedRead& lhs, const AlignedRead& rhs) noexcept;

struct IsDuplicate
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
};

namespace detail {

template <typename Range>
bool no_duplicates(const AlignedRead& read, const Range& reads) noexcept
{
    return std::none_of(std::cbegin(reads), std::cend(reads), [&read] (auto other_itr) { return IsDuplicate{}(read, *other_itr); });
}

struct NextSegmentLess
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
};

} // namespace detail

template <typename ForwardIt>
std::vector<std::vector<ForwardIt>>
find_duplicates(ForwardIt first, const ForwardIt last)
{
    std::vector<std::vector<ForwardIt>> result {};
    std::vector<ForwardIt> buffer {};
    // Recall that reads come sorted w.r.t operator< and it is therefore not guaranteed that 'duplicate'
    // reads (according to IsDuplicate) will be adjacent to one another. In particular, operator< only
    // guarantees that duplicate read segment described in the AlignedRead object will be adjacent.
    const static auto are_primary_dups = [] (const auto& lhs, const auto& rhs) { return primary_segments_are_duplicates(lhs, rhs); };
    for (first = std::adjacent_find(first, last, are_primary_dups); first != last; first = std::adjacent_find(first, last, are_primary_dups)) {
        buffer.reserve(8);
        if (first->has_other_segment()) buffer.push_back(first);
        if (std::next(first)->has_other_segment()) buffer.push_back(std::next(first));
        const AlignedRead& primary_dup_read {*first};
        std::advance(first, 2);
        for (; first != last && primary_segments_are_duplicates(*first, primary_dup_read); ++first) {
            if (first->has_other_segment()) buffer.push_back(first);
        }
        std::sort(buffer.begin(), buffer.end(), [] (ForwardIt lhs, ForwardIt rhs) { return detail::NextSegmentLess{}(*lhs, *rhs); });
        const static auto are_dups = [] (ForwardIt lhs, ForwardIt rhs) { return are_duplicates(*lhs, *rhs); };
        for (auto dup_itr = std::adjacent_find(buffer.cbegin(), buffer.cend(), are_dups); dup_itr != buffer.cend();
                 dup_itr = std::adjacent_find(dup_itr, buffer.cend(), are_dups)) {
            std::vector<ForwardIt> duplicates {dup_itr, std::next(dup_itr, 2)};
            const AlignedRead& duplicate_read {**dup_itr};
            std::advance(dup_itr, 2);
            for (; dup_itr != buffer.cend() && are_duplicates(**dup_itr, duplicate_read); ++dup_itr) {
                duplicates.push_back(*dup_itr);
            }
            result.push_back(std::move(duplicates));
        }
        buffer.clear();
    }
    return result;
}

template <typename ForwardIt>
ForwardIt remove_duplicates(ForwardIt first, ForwardIt last)
{
    // See comment in 'find_duplicates'
    first = std::adjacent_find(first, last, [] (const auto& lhs, const auto& rhs) { return primary_segments_are_duplicates(lhs, rhs); });
    if (first != last) {
        std::vector<ForwardIt> buffer {first++};
        buffer.reserve(100);
        for (auto itr = first; itr != last; ++itr) {
            if (primary_segments_are_duplicates(*itr, *buffer.front())) { // can check any read in buffer
                if (detail::no_duplicates(*itr, buffer)) {
                    if (itr != first) *first = std::move(*itr);
                    buffer.emplace_back(first++);
                }
            } else {
                if (itr != first) *first = std::move(*itr);
                buffer.assign({first++});
            }
        }
    }
    return first;
}

} // namespace

#endif //read_duplicates_hpp
