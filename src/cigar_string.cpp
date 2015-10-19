//
//  cigar_string.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cigar_string.hpp"

#include <ctype.h>   // std::isdigit
#include <algorithm> // std::copy

CigarOperation::CigarOperation(SizeType size, char flag) noexcept
:
size_ {size},
flag_ {flag}
{}

CigarOperation::SizeType CigarOperation::get_size() const noexcept
{
    return size_;
}

char CigarOperation::get_flag() const noexcept
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
    result.reserve(cigar_string.size() / 2);
    std::string digits {};
    
    for (char c : cigar_string) {
        if (std::isdigit(c)) {
            digits += c;
        } else {
            result.emplace_back(std::stoi(digits), c);
            digits.clear();
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

bool is_front_soft_clipped(const CigarString& cigar_string) noexcept
{
    return cigar_string.empty() || cigar_string.front().get_flag() == CigarOperation::SOFT_CLIPPED;
}

bool is_back_soft_clipped(const CigarString& cigar_string) noexcept
{
    return cigar_string.empty() || cigar_string.back().get_flag() == CigarOperation::SOFT_CLIPPED;
}

bool is_soft_clipped(const CigarString& cigar_string) noexcept
{
    return is_front_soft_clipped(cigar_string) || is_back_soft_clipped(cigar_string);
}

std::pair<CigarOperation::SizeType, CigarOperation::SizeType>
get_soft_clipped_sizes(const CigarString& cigar_string) noexcept
{
    if (!is_soft_clipped(cigar_string)) {
        return {0, 0};
    } else {
        auto front_soft_clipped_size = (is_front_soft_clipped(cigar_string)) ?
        cigar_string.front().get_size() : 0;
        auto back_soft_clipped_size = (is_back_soft_clipped(cigar_string)) ?
        cigar_string.back().get_size() : 0;
        return {front_soft_clipped_size, back_soft_clipped_size};
    }
}

// non-member functions

std::ostream& operator<<(std::ostream& os, const CigarOperation& cigar_operation)
{
    os << cigar_operation.get_size() << cigar_operation.get_flag();
    return os;
}

std::ostream& operator<<(std::ostream& os, const CigarString& cigar_string)
{
    std::copy(std::cbegin(cigar_string), std::cend(cigar_string), std::ostream_iterator<CigarOperation>(os));
    return os;
}
