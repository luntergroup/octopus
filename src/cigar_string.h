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
#include <algorithm> // std::copy
#include <numeric>   // std::accumulate
#include <ctype.h>   // std::isdigit
#include <boost/functional/hash.hpp> // boost::hash_combine, boost::hash_range

#include "equitable.h"

class CigarOperation : public Comparable<CigarOperation> // Comparable so can compare reads
{
public:
    using SizeType = std::uint_fast32_t;
    
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
    
    CigarOperation() = delete;
    explicit CigarOperation(SizeType size, char type) noexcept;
    
    CigarOperation(const CigarOperation&)            = default;
    CigarOperation& operator=(const CigarOperation&) = default;
    CigarOperation(CigarOperation&&)                 = default;
    CigarOperation& operator=(CigarOperation&&)      = default;
    
    SizeType get_size() const noexcept;
    char get_flag() const noexcept;
    bool advances_reference() const noexcept;
    bool advances_sequence() const noexcept;
    
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

inline bool CigarOperation::advances_reference() const noexcept
{
    return !(flag_ == INSERTION || flag_ == HARD_CLIPPED || flag_ == PADDING);
}

inline bool CigarOperation::advances_sequence() const noexcept
{
    return !(flag_ == DELETION || flag_ == HARD_CLIPPED);
}

inline CigarString parse_cigar_string(const std::string& cigar_string)
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

inline bool is_front_soft_clipped(const CigarString& cigar_string) noexcept
{
    return (cigar_string.size() > 0) ? cigar_string.front().get_flag() == 'S' : false;
}

inline bool is_back_soft_clipped(const CigarString& cigar_string) noexcept
{
    return (cigar_string.size() > 0) ? cigar_string.back().get_flag() == 'S' : false;
}

inline bool is_soft_clipped(const CigarString& cigar_string) noexcept
{
    return is_front_soft_clipped(cigar_string) || is_back_soft_clipped(cigar_string);
}

inline std::pair<CigarOperation::SizeType, CigarOperation::SizeType>
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

template <typename SizeType>
inline
SizeType soft_clipped_read_begin(const CigarString& cigar_string, SizeType hard_clipped_begin) noexcept
{
    if (is_front_soft_clipped(cigar_string)) {
        hard_clipped_begin -= static_cast<SizeType>(cigar_string.at(0).get_size());
    }
    return hard_clipped_begin;
}

template <typename SizeType=unsigned>
inline SizeType operations_size(const CigarString& cigar_string) noexcept
{
    return std::accumulate(std::cbegin(cigar_string), std::cend(cigar_string), SizeType {},
                           [] (const unsigned lhs, const CigarOperation& rhs) {
                               return lhs + rhs.get_size();
                           });
}

template <typename SizeType=unsigned>
inline SizeType reference_size(const CigarString& cigar_string) noexcept
{
    SizeType result {};
    
    for (const auto& op : cigar_string) {
        if (op.advances_reference()) {
            result += op.get_size();
        }
    }
    
    return result;
}

template <typename SizeType=unsigned>
inline SizeType sequence_size(const CigarString& cigar_string) noexcept
{
    SizeType result {};
    
    for (const auto& op : cigar_string) {
        if (op.advances_sequence()) {
            result += op.get_size();
        }
    }
    
    return result;
}

template <typename SizeType>
CigarOperation get_operation_at_sequence_index(const CigarString& cigar_string, SizeType the_sequence_index)
{
    auto first = std::cbegin(cigar_string);
    
    while (the_sequence_index >= first->get_size()) {
        ++first;
        the_sequence_index -= first->get_size();
    }
    
    return *first;
}

namespace detail
{
    inline bool advances_reference(const CigarOperation& op)
    {
        return op.advances_reference();
    }
    
    inline bool advances_sequence(const CigarOperation& op)
    {
        return op.advances_sequence();
    }
    
    template <typename SizeType, typename Predicate>
    CigarString splice(const CigarString& cigar_string, SizeType offset, SizeType size, Predicate pred)
    {
        CigarString result {};
        result.reserve(cigar_string.size()); // ensures no reallocation due to emplace_back
        
        auto first = std::cbegin(cigar_string);
        auto last  = std::cend(cigar_string);
        
        while (first != last && (offset >= first->get_size() || !pred(*first))) {
            if (pred(*first)) offset -= first->get_size();
            ++first;
        }
        
        if (first != last) {
            auto remainder = first->get_size() - offset;
            
            if (remainder >= size) {
                result.emplace_back(size, first->get_flag());
                result.shrink_to_fit();
                return result;
            }
            
            result.emplace_back(remainder, first->get_flag());
            size -= remainder;
            ++first;
        }
        
        while (first != last && (size >= first->get_size() || !pred(*first))) {
            result.emplace_back(*first);
            if (pred(*first)) size -= first->get_size();
            ++first;
        }
        
        if (first != last && size > 0) {
            result.emplace_back(size, first->get_flag());
        }
        
        result.shrink_to_fit();
        
        return result;
    }
} // end namespace detail

template <typename SizeType>
CigarString operation_splice(const CigarString& cigar_string, SizeType offset, SizeType size)
{
    return detail::splice(cigar_string, offset, size, [] (const auto& op) { return true; });
}

template <typename SizeType>
CigarString reference_splice(const CigarString& cigar_string, SizeType offset, SizeType size)
{
    return detail::splice(cigar_string, offset, size, [] (const auto& op) { return detail::advances_reference(op); });
}

template <typename SizeType>
CigarString sequence_splice(const CigarString& cigar_string, SizeType offset, SizeType size)
{
    return detail::splice(cigar_string, offset, size, [] (const auto& op) { return detail::advances_sequence(op); });
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

// Note CigarString gets all of std::vector comparison methods. In particular, operator< uses
// std::lexicographical_compare

inline std::ostream& operator<<(std::ostream& os, const CigarOperation& cigar_operation)
{
    os << cigar_operation.get_size() << cigar_operation.get_flag();
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const CigarString& cigar_string)
{
    std::copy(std::cbegin(cigar_string), std::cend(cigar_string), std::ostream_iterator<CigarOperation>(os));
    return os;
}

namespace std {
    template <> struct hash<CigarOperation>
    {
        size_t operator()(const CigarOperation& op) const
        {
            size_t seed {};
            boost::hash_combine(seed, op.get_flag());
            boost::hash_combine(seed, op.get_size());
            return seed;
        }
    };
}

namespace boost {
    template <> struct hash<CigarOperation> : std::unary_function<CigarOperation, std::size_t>
    {
        std::size_t operator()(const CigarOperation& op) const
        {
            return std::hash<CigarOperation>()(op);
        }
    };
}

namespace std {
    template <> struct hash<CigarString>
    {
        size_t operator()(const CigarString& cigar) const
        {
            return boost::hash_range(std::cbegin(cigar), std::cend(cigar));
        }
    };
}

#endif
