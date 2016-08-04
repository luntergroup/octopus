// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cigar_string.hpp"

#include <cctype>
#include <algorithm>
#include <array>
#include <ostream>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

namespace octopus {

constexpr char CigarOperation::ALIGNMENT_MATCH;
constexpr char CigarOperation::SEQUENCE_MATCH;
constexpr char CigarOperation::SUBSTITUTION;
constexpr char CigarOperation::INSERTION;
constexpr char CigarOperation::DELETION;
constexpr char CigarOperation::SOFT_CLIPPED;
constexpr char CigarOperation::HARD_CLIPPED;
constexpr char CigarOperation::PADDING;
constexpr char CigarOperation::SKIPPED;

CigarOperation::CigarOperation(const Size size, const Flag flag) noexcept
:
size_ {size},
flag_ {flag}
{}

CigarOperation::Size CigarOperation::size() const noexcept
{
    return size_;
}

CigarOperation::Flag CigarOperation::flag() const noexcept
{
    return flag_;
}

bool CigarOperation::advances_reference() const noexcept
{
    return !(flag_ == INSERTION || flag_ == HARD_CLIPPED || flag_ == PADDING);
}

bool CigarOperation::advances_sequence() const noexcept
{
    return !(flag_ == DELETION || flag_ == HARD_CLIPPED);
}

bool is_match(const CigarOperation& op) noexcept
{
    switch (op.flag()) {
        case CigarOperation::ALIGNMENT_MATCH:
        case CigarOperation::SEQUENCE_MATCH:
        case CigarOperation::SUBSTITUTION: return true;
        default: return false;
    }
}

bool is_indel(const CigarOperation& op) noexcept
{
    return op.flag() == CigarOperation::INSERTION || op.flag() == CigarOperation::DELETION;
}

bool is_clipping(const CigarOperation& op) noexcept
{
    return op.flag() == CigarOperation::SOFT_CLIPPED || op.flag() == CigarOperation::HARD_CLIPPED;
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
            result.emplace_back(boost::lexical_cast<CigarOperation::Size>(digits), c);
            digits.clear();
        }
    }
    
    if (!digits.empty()) {
        throw std::invalid_argument {"parse_cigar: unparsed characters in " + cigar};
    }
    
    result.shrink_to_fit();
    
    return result;
}

bool is_valid_flag(const CigarOperation& op) noexcept
{
    static constexpr std::array<char, 9> valid_flags {
        CigarOperation::ALIGNMENT_MATCH,
        CigarOperation::INSERTION,
        CigarOperation::DELETION,
        CigarOperation::SOFT_CLIPPED,
        CigarOperation::HARD_CLIPPED,
        CigarOperation::SKIPPED,
        CigarOperation::SEQUENCE_MATCH,
        CigarOperation::SUBSTITUTION,
        CigarOperation::PADDING
    };
    
    return std::find(std::cbegin(valid_flags), std::cend(valid_flags), op.flag()) != std::cend(valid_flags);
}

bool is_valid(const CigarString& cigar) noexcept
{
    return !cigar.empty() && std::all_of(std::cbegin(cigar), std::cend(cigar),
                                         [] (const auto& op) {
                                             return op.size() > 0 && is_valid_flag(op);
                                         });
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
    return !cigar.empty() && cigar.front().flag() == CigarOperation::SOFT_CLIPPED;
}

bool is_back_soft_clipped(const CigarString& cigar) noexcept
{
    return !cigar.empty() && cigar.back().flag() == CigarOperation::SOFT_CLIPPED;
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

template <typename Predicate>
CigarString splice(const CigarString& cigar, CigarOperation::Size offset,
                   CigarOperation::Size size, Predicate pred)
{
    CigarString result {};
    result.reserve(cigar.size());
    
    auto op_it = std::cbegin(cigar);
    
    const auto last_op = std::cend(cigar);
    
    while (op_it != last_op && (offset >= op_it->size() || !pred(*op_it))) {
        if (pred(*op_it)) offset -= op_it->size();
        ++op_it;
    }
    
    if (op_it != last_op) {
        const auto remainder = op_it->size() - offset;
        
        if (remainder >= size) {
            result.emplace_back(size, op_it->flag());
            result.shrink_to_fit();
            return result;
        }
        
        result.emplace_back(remainder, op_it->flag());
        size -= remainder;
        ++op_it;
    }
    
    while (op_it != last_op && size > 0 && (size >= op_it->size() || !pred(*op_it))) {
        result.emplace_back(*op_it);
        if (pred(*op_it)) size -= op_it->size();
        ++op_it;
    }
    
    if (op_it != last_op && size > 0) {
        result.emplace_back(size, op_it->flag());
    }
    
    result.shrink_to_fit();
    
    return result;
}

CigarString splice(const CigarString& cigar, const CigarOperation::Size offset,
                   const CigarOperation::Size size)
{
    return splice(cigar, offset, size, [] (const auto& op) { return true; });
}

CigarString splice_reference(const CigarString& cigar, const CigarOperation::Size offset,
                             const CigarOperation::Size size)
{
    return splice(cigar, offset, size, [] (const auto& op) { return op.advances_reference(); });
}

CigarString splice_sequence(const CigarString& cigar, const CigarOperation::Size offset,
                            const CigarOperation::Size size)
{
    return splice(cigar, offset, size, [] (const auto& op) { return op.advances_sequence(); });
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
    using boost::hash_combine;
    size_t result {};
    hash_combine(result, op.flag());
    hash_combine(result, op.size());
    return result;
}

std::size_t CigarHash::operator()(const CigarString& cigar) const noexcept
{
    return boost::hash_range(std::cbegin(cigar), std::cend(cigar));
}

} // namespace octopus
