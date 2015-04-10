//
//  cigar_string.h
//  Octopus
//
//  Created by Daniel Cooke on 12/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_cigar_string_h
#define Octopus_cigar_string_h

#include <string>
#include <cstdint>
#include <iterator>
#include <vector>
#include <ostream>
#include <algorithm> // std::copy, std::minmax, std::mismatch
#include <ctype.h>   // std::isdigit

#include "equitable.h"

class CigarOperation : public Comparable<CigarOperation> // Comparable so can compare reads
{
public:
    using SizeType = std::uint_fast32_t;
    
    CigarOperation() = delete;
    explicit CigarOperation(SizeType size, char type) noexcept;
    
    CigarOperation(const CigarOperation&)            = default;
    CigarOperation& operator=(const CigarOperation&) = default;
    CigarOperation(CigarOperation&&)                 = default;
    CigarOperation& operator=(CigarOperation&&)      = default;
    
    SizeType get_size() const noexcept;
    char get_flag() const noexcept;
    
    // No getting around these really - they have to go somewhere. Could go into a config,
    // but would loose constexpr-ness.
    static const constexpr char ALIGNMENT_MATCH {'M'};
    static const constexpr char SEQUENCE_MATCH  {'='};
    static const constexpr char SUBSTITUTION    {'X'};
    static const constexpr char INSERTION       {'I'};
    static const constexpr char DELETION        {'D'};
    static const constexpr char SOFT_CLIPPED    {'S'};
    static const constexpr char HARD_CLIPPED    {'H'};
    static const constexpr char PADDING         {'P'};
    static const constexpr char SKIPPED         {'N'};
    
private:
    SizeType size_;
    char flag_;
};

using CigarString = std::vector<CigarOperation>;

inline CigarOperation::CigarOperation(SizeType size, char flag) noexcept
:
size_ {size},
flag_ {flag}
{}

inline CigarOperation::SizeType CigarOperation::get_size() const noexcept
{
    return size_;
}

inline char CigarOperation::get_flag() const noexcept
{
    return flag_;
}

inline CigarString parse_cigar_string(const std::string& a_cigar_string)
{
    CigarString result {};
    result.reserve(a_cigar_string.size() / 2);
    std::string digits {};
    
    for (char c : a_cigar_string) {
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

inline bool is_front_soft_clipped(const CigarString& a_cigar_string) noexcept
{
    return (a_cigar_string.size() > 0) ? a_cigar_string.front().get_flag() == 'S' : false;
}

inline bool is_back_soft_clipped(const CigarString& a_cigar_string) noexcept
{
    return (a_cigar_string.size() > 0) ? a_cigar_string.back().get_flag() == 'S' : false;
}

inline bool is_soft_clipped(const CigarString& a_cigar_string) noexcept
{
    return is_front_soft_clipped(a_cigar_string) || is_back_soft_clipped(a_cigar_string);
}

inline std::pair<CigarOperation::SizeType, CigarOperation::SizeType>
get_soft_clipped_sizes(const CigarString& a_cigar_string) noexcept
{
    if (!is_soft_clipped(a_cigar_string)) {
        return {0, 0};
    } else {
        auto front_soft_clipped_size = (is_front_soft_clipped(a_cigar_string)) ?
            a_cigar_string.front().get_size() : 0;
        auto back_soft_clipped_size = (is_back_soft_clipped(a_cigar_string)) ?
            a_cigar_string.back().get_size() : 0;
        return {front_soft_clipped_size, back_soft_clipped_size};
    }
}

template <typename T>
inline T get_soft_clipped_read_begin(const CigarString& a_cigar_string, T hard_clipped_begin) noexcept
{
    if (is_front_soft_clipped(a_cigar_string)) {
        hard_clipped_begin -= static_cast<T>(a_cigar_string.at(0).get_size());
    }
    return hard_clipped_begin;
}

inline bool operator==(const CigarOperation& lhs, const CigarOperation& rhs)
{
    return lhs.get_flag() == rhs.get_flag() && lhs.get_size() == rhs.get_size();
}

inline bool operator<(const CigarOperation& lhs, const CigarOperation& rhs)
{
    return (lhs.get_flag() == rhs.get_flag()) ? lhs.get_size() < rhs.get_size() :
                                                lhs.get_flag() < rhs.get_flag();
}

inline bool operator==(const CigarString& lhs, const CigarString& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}
inline bool operator!=(const CigarString& lhs, const CigarString& rhs){return !operator==(lhs, rhs);}
inline bool operator<(const CigarString& lhs, const CigarString& rhs)
{
    auto p = std::minmax(lhs, rhs, [] (const auto& lhs, const auto& rhs) { return lhs.size() <= rhs.size(); });
    auto itrs = std::mismatch(std::cbegin(p.first), std::cend(p.first), std::cbegin(p.second));
    return (itrs.first != std::cend(p.first)) ? *(itrs.first) < *(itrs.second) : false;
}
inline bool operator> (const CigarString& lhs, const CigarString& rhs){return  operator< (rhs,lhs);}
inline bool operator<=(const CigarString& lhs, const CigarString& rhs){return !operator> (lhs,rhs);}
inline bool operator>=(const CigarString& lhs, const CigarString& rhs){return !operator< (lhs,rhs);}

inline std::ostream& operator<<(std::ostream& os, const CigarOperation& a_cigar_operation)
{
    os << a_cigar_operation.get_size() << a_cigar_operation.get_flag();
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const CigarString& a_cigar_string)
{
    std::copy(std::cbegin(a_cigar_string), std::cend(a_cigar_string),
              std::ostream_iterator<CigarOperation>(os));
    return os;
}

#endif
