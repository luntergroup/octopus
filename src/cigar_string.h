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

#include "equitable.h"

class CigarOperation : Equitable<CigarOperation>
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
    
    // No getting arounf these really - they have to go somewhere. Could go into a config,
    // but would loose constexpr-ness.
    static const constexpr char ALIGNMENT_MATCH {'M'};
    static const constexpr char SEQUENCE_MATCH {'='};
    static const constexpr char MISMATCH {'X'};
    static const constexpr char INSERTION {'I'};
    static const constexpr char DELETION {'D'};
    static const constexpr char SOFT_CLIPPED {'S'};
    static const constexpr char HARD_CLIPPED {'H'};
    static const constexpr char PADDING {'P'};
    static const constexpr char SKIPPED {'N'};
    
private:
    SizeType size_;
    char flag_;
};

using CigarString = std::vector<CigarOperation>;

inline CigarOperation::CigarOperation(SizeType size, char flag) noexcept
:size_ {size},
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
    const auto num_operations = a_cigar_string.size() / 2;
    CigarString result {};
    result.reserve(num_operations);
    for (unsigned i {0}; i < num_operations; ++i) {
        result.emplace_back(a_cigar_string[i], a_cigar_string[i + 1]);
    }
    return result;
}

inline bool is_soft_clipped(const CigarString& a_cigar_string) noexcept
{
    return (a_cigar_string.size() > 0) ? a_cigar_string[0].get_flag() == 'S' : false;
}

template <typename T>
inline T get_soft_clipped_read_begin(const CigarString& a_cigar_string, T hard_clipped_begin) noexcept
{
    if (is_soft_clipped(a_cigar_string)) {
        hard_clipped_begin -= static_cast<T>(a_cigar_string.at(0).get_size());
    }
    return hard_clipped_begin;
}

inline bool operator==(const CigarOperation& lhs, const CigarOperation& rhs)
{
    return lhs.get_size() == rhs.get_size() && lhs.get_size() == rhs.get_size();
}

inline bool operator==(const CigarString& lhs, const CigarString& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

inline bool operator!=(const CigarString& lhs, const CigarString& rhs)
{
    return !operator==(lhs, rhs);
}

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
