// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cigar_string.hpp"

#include <cctype>
#include <algorithm>
#include <array>
#include <ostream>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include "utils/append.hpp"

namespace octopus {

CigarOperation::CigarOperation(const Size size, const Flag flag) noexcept
: size_ {size}
, flag_ {flag}
{}

void CigarOperation::set_flag(Flag type) noexcept
{
    flag_ = type;
}

void CigarOperation::set_size(Size size) noexcept
{
    size_ = size;
}

CigarOperation::Flag CigarOperation::flag() const noexcept
{
    return flag_;
}

CigarOperation::Size CigarOperation::size() const noexcept
{
    return size_;
}

// non-member methods

bool is_valid(const CigarOperation::Flag flag) noexcept
{
    using Flag = CigarOperation::Flag;
    switch (flag) {
        case Flag::alignmentMatch:
        case Flag::sequenceMatch:
        case Flag::substitution:
        case Flag::insertion:
        case Flag::deletion:
        case Flag::softClipped:
        case Flag::hardClipped:
        case Flag::padding:
        case Flag::skipped: return true;
        default: return false;
    }
}

bool is_valid(const CigarOperation& op) noexcept
{
    return is_valid(op.flag()) && op.size() > 0;
}

bool advances_reference(CigarOperation::Flag flag) noexcept
{
    using Flag = CigarOperation::Flag;
    return !(flag == Flag::insertion || flag == Flag::hardClipped || flag == Flag::padding);
}

bool advances_reference(const CigarOperation& op) noexcept
{
    return advances_reference(op.flag());
}

bool advances_sequence(CigarOperation::Flag flag) noexcept
{
    using Flag = CigarOperation::Flag;
    return !(flag == Flag::deletion || flag == Flag::hardClipped);
}

bool advances_sequence(const CigarOperation& op) noexcept
{
    return advances_sequence(op.flag());
}

bool is_match(CigarOperation::Flag flag) noexcept
{
    using Flag = CigarOperation::Flag;
    switch (flag) {
        case Flag::alignmentMatch:
        case Flag::sequenceMatch:
        case Flag::substitution: return true;
        default: return false;
    }
}

bool is_match(const CigarOperation& op) noexcept
{
    return is_match(op.flag());
}

bool is_insertion(CigarOperation::Flag flag) noexcept
{
    return flag == CigarOperation::Flag::insertion;
}

bool is_insertion(const CigarOperation& op) noexcept
{
    return is_insertion(op.flag());
}

bool is_deletion(CigarOperation::Flag flag) noexcept
{
    return flag == CigarOperation::Flag::deletion;
}

bool is_deletion(const CigarOperation& op) noexcept
{
    return is_deletion(op.flag());
}

bool is_indel(CigarOperation::Flag flag) noexcept
{
    return is_insertion(flag) || is_deletion(flag);
}

bool is_indel(const CigarOperation& op) noexcept
{
    return is_indel(op.flag());
}

bool is_clipping(CigarOperation::Flag flag) noexcept
{
    using Flag = CigarOperation::Flag;
    return flag == Flag::softClipped || flag == Flag::hardClipped;
}

bool is_clipping(const CigarOperation& op) noexcept
{
    return is_clipping(op.flag());
}

// CigarString

CigarString parse_cigar(const std::string& cigar)
{
    CigarString result {};
    result.reserve(cigar.size() / 2); // max possible CigarOperation
    std::string digits {};
    digits.reserve(3); // 100+ matches are common
    
    for (const char c : cigar) {
        if (std::isdigit(c)) {
            digits += c;
        } else {
            result.emplace_back(boost::lexical_cast<CigarOperation::Size>(digits),
                                static_cast<CigarOperation::Flag>(c));
            digits.clear();
        }
    }
    
    if (!digits.empty()) {
        throw std::invalid_argument {"parse_cigar: unparsed characters in " + cigar};
    }
    
    result.shrink_to_fit();
    
    return result;
}

bool is_valid(const CigarString& cigar) noexcept
{
    return !cigar.empty() && std::all_of(std::cbegin(cigar), std::cend(cigar),
                                         [] (const auto& op) { return is_valid(op); });
}

bool is_minimal(const CigarString& cigar) noexcept
{
    return std::adjacent_find(std::cbegin(cigar), std::cend(cigar),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.flag() == rhs.flag();
                              }) == std::cend(cigar);
}

bool is_front_soft_clipped(const CigarString& cigar) noexcept
{
    return !cigar.empty() && cigar.front().flag() == CigarOperation::Flag::softClipped;
}

bool is_back_soft_clipped(const CigarString& cigar) noexcept
{
    return !cigar.empty() && cigar.back().flag() == CigarOperation::Flag::softClipped;
}

bool is_soft_clipped(const CigarString& cigar) noexcept
{
    return is_front_soft_clipped(cigar) || is_back_soft_clipped(cigar);
}

std::pair<CigarOperation::Size, CigarOperation::Size>
get_soft_clipped_sizes(const CigarString& cigar) noexcept
{
    return std::make_pair((is_front_soft_clipped(cigar)) ? cigar.front().size() : 0,
                          (is_back_soft_clipped(cigar)) ? cigar.back().size() : 0);
}

// non-member functions

namespace {

template <typename Predicate1, typename Predicate2>
CigarString copy(const CigarString& cigar, CigarOperation::Size offset, CigarOperation::Size size,
                 Predicate1 offset_pred, Predicate2 size_pred)
{
    CigarString result {};
    result.reserve(cigar.size());
    auto op_itr = std::cbegin(cigar);
    const auto last_op_itr = std::cend(cigar);
    while (op_itr != last_op_itr && offset > 0 && (offset >= op_itr->size() || !offset_pred(*op_itr))) {
        if (offset_pred(*op_itr)) {
            offset -= op_itr->size();
        }
        ++op_itr;
    }
    if (op_itr != last_op_itr) {
        const auto remainder = op_itr->size() - offset;
        if (remainder >= size) {
            if (size_pred(*op_itr)) {
                result.emplace_back(size, op_itr->flag());
            } else {
                result.emplace_back(op_itr->size(), op_itr->flag());
            }
            return result;
        }
        result.emplace_back(remainder, op_itr->flag());
        size -= remainder;
        ++op_itr;
    }
    while (op_itr != last_op_itr && size > 0 && (size >= op_itr->size() || !size_pred(*op_itr))) {
        result.emplace_back(*op_itr);
        if (size_pred(*op_itr)) {
            size -= op_itr->size();
        }
        ++op_itr;
    }
    if (op_itr != last_op_itr && size > 0) {
        result.emplace_back(size, op_itr->flag());
    }
    return result;
}

template <typename Predicate>
CigarString copy(const CigarString& cigar, CigarOperation::Size offset, CigarOperation::Size size,
                 Predicate pred)
{
    return copy(cigar, offset, size, pred, pred);
}

struct AdvancesReferencePred
{
    bool operator()(const CigarOperation& op) const noexcept
    {
        return advances_reference(op);
    }
};
struct AdvancesSequencePred
{
    bool operator()(const CigarOperation& op) const noexcept
    {
        return advances_sequence(op);
    }
};

template <typename Predicate>
CigarString copy(const CigarString& cigar, CigarOperation::Size offset, CigarOperation::Size size,
                 Predicate offset_pred, const CigarStringCopyPolicy size_policy)
{
    using CopyPolicy = CigarStringCopyPolicy;
    switch (size_policy) {
        case CopyPolicy::reference:
            return copy(cigar, offset, size, offset_pred, AdvancesReferencePred {});
        case CopyPolicy::sequence:
            return copy(cigar, offset, size, offset_pred, AdvancesSequencePred {});
        case CopyPolicy::both:
        default:
            return copy(cigar, offset, size, offset_pred, [] (const auto& op) { return true; });
    }
}

} // namespace

CigarString copy(const CigarString& cigar, CigarOperation::Size offset, CigarOperation::Size size,
                 const CigarStringCopyPolicy offset_policy, const CigarStringCopyPolicy size_policy)
{
    using CopyPolicy = CigarStringCopyPolicy;
    switch (offset_policy) {
        case CopyPolicy::reference:
            return copy(cigar, offset, size, AdvancesReferencePred {}, size_policy);
        case CopyPolicy::sequence:
            return copy(cigar, offset, size, AdvancesSequencePred {}, size_policy);
        case CopyPolicy::both:
        default:
            return copy(cigar, offset, size, [] (const auto& op) { return true; }, size_policy);
    }
}

CigarString copy_reference(const CigarString& cigar, CigarOperation::Size offset, CigarOperation::Size size)
{
    return copy(cigar, offset, size, AdvancesReferencePred {});
}

CigarString copy_sequence(const CigarString& cigar, CigarOperation::Size offset, CigarOperation::Size size)
{
    return copy(cigar, offset, size, AdvancesSequencePred {});
}

std::vector<CigarOperation::Flag> decompose(const CigarString& cigar)
{
    std::vector<CigarOperation::Flag> result {};
    result.reserve(sum_operation_sizes(cigar));
    for (const auto& op : cigar) {
        utils::append(result, op.size(), op.flag());
    }
    return result;
}

CigarString collapse_matches(const CigarString& cigar)
{
    CigarString result {};
    result.reserve(cigar.size());
    for (auto match_end_itr = std::begin(cigar); match_end_itr != std::cend(cigar); ) {
        const auto f_is_match = [] (const CigarOperation& op) { return is_match(op); };
        const auto match_begin_itr = std::find_if(match_end_itr, std::end(cigar), f_is_match);
        result.insert(std::cend(result), match_end_itr, match_begin_itr);
        if (match_begin_itr == std::cend(cigar)) break;
        match_end_itr = std::find_if_not(std::next(match_begin_itr), std::end(cigar), f_is_match);
        auto match_size = std::accumulate(match_begin_itr, match_end_itr, 0,
                                          [] (auto curr, const auto& op) { return curr + op.size(); });
        result.emplace_back(match_size, CigarOperation::Flag::alignmentMatch);
    }
    return result;
}

bool operator==(const CigarOperation& lhs, const CigarOperation& rhs) noexcept
{
    return lhs.flag() == rhs.flag() && lhs.size() == rhs.size();
}

bool operator<(const CigarOperation& lhs, const CigarOperation& rhs) noexcept
{
    return (lhs.flag() == rhs.flag()) ? lhs.size() < rhs.size() : lhs.flag() < rhs.flag();
}

std::ostream& operator<<(std::ostream& os, const CigarOperation::Flag& flag)
{
    os << static_cast<char>(flag);
    return os;
}

std::ostream& operator<<(std::ostream& os, const CigarOperation& cigar_operation)
{
    os << cigar_operation.size() << cigar_operation.flag();
    return os;
}

std::ostream& operator<<(std::ostream& os, const CigarString& cigar)
{
    std::copy(std::cbegin(cigar), std::cend(cigar), std::ostream_iterator<CigarOperation>(os));
    return os;
}

std::size_t CigarHash::operator()(const CigarOperation& op) const noexcept
{
    std::size_t result {};
    using boost::hash_combine;
    hash_combine(result, op.flag());
    hash_combine(result, op.size());
    return result;
}

std::size_t CigarHash::operator()(const CigarString& cigar) const noexcept
{
    return boost::hash_range(std::cbegin(cigar), std::cend(cigar));
}

} // namespace octopus
