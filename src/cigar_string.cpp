//
//  cigar_string.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cigar_string.hpp"

#include <cctype>
#include <algorithm>
#include <array>
#include <ostream>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

CigarOperation::CigarOperation(const SizeType size, const char flag) noexcept
:
size_ {size},
flag_ {flag}
{}

CigarOperation::SizeType CigarOperation::size() const noexcept
{
    return size_;
}

char CigarOperation::flag() const noexcept
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

CigarString parse_cigar_string(const std::string& cigar_string)
{
    CigarString result {};
    result.reserve(cigar_string.size() / 2); // max possible CigarOperation
    
    std::string digits {};
    digits.reserve(3); // 100+ matches are common
    
    for (const char c : cigar_string) {
        if (std::isdigit(c)) {
            digits += c;
        } else {
            try {
                result.emplace_back(boost::lexical_cast<CigarOperation::SizeType>(digits), c);
                digits.clear();
            } catch (const boost::bad_lexical_cast&) {
                throw;
            }
        }
    }
    
    if (!digits.empty()) {
        throw std::invalid_argument {"parse_cigar_string: could not parse all characters of " + cigar_string};
    }
    
    result.shrink_to_fit();
    
    return result;
}

bool is_valid_flag(const CigarOperation& op)
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
    
    static const auto first_valid = std::cbegin(valid_flags);
    static const auto last_valid  = std::cend(valid_flags);
    
    return std::find(first_valid, last_valid, op.flag()) != last_valid;
}

bool is_valid_cigar(const CigarString& cigar) noexcept
{
    return !cigar.empty() && std::all_of(std::cbegin(cigar), std::cend(cigar),
                                         [] (const auto& op) {
                                             return op.size() > 0 && is_valid_flag(op);
                                         });
}

bool is_minimal_cigar(const CigarString& cigar) noexcept
{
    return std::adjacent_find(std::cbegin(cigar), std::cend(cigar),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.flag() == rhs.flag();
                              }) == std::cend(cigar);
}

bool is_front_soft_clipped(const CigarString& cigar_string) noexcept
{
    return cigar_string.empty() || cigar_string.front().flag() == CigarOperation::SOFT_CLIPPED;
}

bool is_back_soft_clipped(const CigarString& cigar_string) noexcept
{
    return cigar_string.empty() || cigar_string.back().flag() == CigarOperation::SOFT_CLIPPED;
}

bool is_soft_clipped(const CigarString& cigar_string) noexcept
{
    return is_front_soft_clipped(cigar_string) || is_back_soft_clipped(cigar_string);
}

std::pair<CigarOperation::SizeType, CigarOperation::SizeType>
get_soft_clipped_sizes(const CigarString& cigar_string) noexcept
{
    return std::make_pair((is_front_soft_clipped(cigar_string)) ? cigar_string.front().size() : 0,
                          (is_back_soft_clipped(cigar_string)) ? cigar_string.back().size() : 0);
}

// non-member functions

template <typename Predicate>
CigarString splice(const CigarString& cigar_string, CigarOperation::SizeType offset,
                   CigarOperation::SizeType size, Predicate pred)
{
    CigarString result {};
    result.reserve(cigar_string.size()); // ensures no reallocation due to emplace_back
    
    auto op_it = std::cbegin(cigar_string);
    
    const auto last = std::cend(cigar_string);
    
    while (op_it != last && (offset >= op_it->size() || !pred(*op_it))) {
        if (pred(*op_it)) offset -= op_it->size();
        ++op_it;
    }
    
    if (op_it != last) {
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
    
    while (op_it != last && size > 0 && (size >= op_it->size() || !pred(*op_it))) {
        result.emplace_back(*op_it);
        if (pred(*op_it)) size -= op_it->size();
        ++op_it;
    }
    
    if (op_it != last && size > 0) {
        result.emplace_back(size, op_it->flag());
    }
    
    result.shrink_to_fit();
    
    return result;
}

CigarString splice(const CigarString& cigar_string, const CigarOperation::SizeType offset,
                   const CigarOperation::SizeType size)
{
    return splice(cigar_string, offset, size, [] (const auto& op) { return true; });
}

CigarString splice(const CigarString& cigar_string, const CigarOperation::SizeType size)
{
    return splice(cigar_string, 0, size);
}

CigarString splice_reference(const CigarString& cigar_string, const CigarOperation::SizeType offset,
                             const CigarOperation::SizeType size)
{
    return splice(cigar_string, offset, size, [] (const auto& op) { return op.advances_reference(); });
}

CigarString splice_reference(const CigarString& cigar_string, const CigarOperation::SizeType size)
{
    return splice(cigar_string, 0, size);
}

CigarString splice_sequence(const CigarString& cigar_string, const CigarOperation::SizeType offset,
                            const CigarOperation::SizeType size)
{
    return splice(cigar_string, offset, size, [] (const auto& op) { return op.advances_sequence(); });
}

CigarString splice_sequence(const CigarString& cigar_string, const CigarOperation::SizeType size)
{
    return splice(cigar_string, 0, size);
}

std::ostream& operator<<(std::ostream& os, const CigarOperation& cigar_operation)
{
    os << cigar_operation.size() << cigar_operation.flag();
    return os;
}

std::ostream& operator<<(std::ostream& os, const CigarString& cigar_string)
{
    std::copy(std::cbegin(cigar_string), std::cend(cigar_string), std::ostream_iterator<CigarOperation>(os));
    return os;
}
