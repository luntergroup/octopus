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
#include <algorithm>

#include "equitable.h"

using std::uint_fast32_t;

class CigarString : Equitable<CigarString>
{
public:
    class CigarOperation
    {
    public:
        CigarOperation() = delete;
        explicit CigarOperation(uint_fast32_t size, char type) noexcept;
        
        CigarOperation(const CigarOperation&)            = default;
        CigarOperation& operator=(const CigarOperation&) = default;
        CigarOperation(CigarOperation&&)                 = default;
        CigarOperation& operator=(CigarOperation&&)      = default;
        
        uint_fast32_t get_size() const noexcept;
        char get_flag() const noexcept;
    private:
        uint_fast32_t size_;
        char flag_;
    };
    
    using Iterator      = std::vector<CigarOperation>::iterator;
    using ConstIterator = std::vector<CigarOperation>::const_iterator;
    
    CigarString() = delete;
    CigarString(const std::string& a_cigar_string);
    explicit CigarString(std::vector<CigarOperation>&& a_cigar_string);
    
    CigarString(const CigarString&)            = default;
    CigarString& operator=(const CigarString&) = default;
    CigarString(CigarString&&)                 = default;
    CigarString& operator=(CigarString&&)      = default;
    
    Iterator begin() noexcept;
    Iterator end() noexcept;
    ConstIterator begin() const noexcept;
    ConstIterator end() const noexcept;
    ConstIterator cbegin() const noexcept;
    ConstIterator cend() const noexcept;
    
    CigarOperation get_cigar_operation(uint_fast32_t index) const;
    const CigarOperation& at(uint_fast32_t n) const;
    bool empty() const noexcept;
    uint_fast32_t size() const noexcept;
    bool is_soft_clipped() const noexcept;
    
private:
    std::vector<CigarOperation> the_cigar_string_;
};

inline CigarString::CigarString(const std::string& a_cigar_string) : the_cigar_string_ {}
{
    const uint_fast32_t num_operations {static_cast<uint_fast32_t>(a_cigar_string.length() / 2)};
    the_cigar_string_.reserve(num_operations);
    for (uint_fast32_t i {0}; i < num_operations; ++i) {
        the_cigar_string_.emplace_back(a_cigar_string[i], a_cigar_string[i + 1]);
    }
}

inline CigarString::CigarOperation::CigarOperation(uint_fast32_t size, char flag) noexcept
:size_ {size},
 flag_ {flag}
{}

inline uint_fast32_t CigarString::CigarOperation::get_size() const noexcept
{
    return size_;
}

inline char CigarString::CigarOperation::get_flag() const noexcept
{
    return flag_;
}

inline CigarString::CigarString(std::vector<CigarOperation>&& a_cigar_string)
:the_cigar_string_ {std::move(a_cigar_string)}
{}

inline CigarString::CigarOperation CigarString::get_cigar_operation(uint_fast32_t index) const
{
    return the_cigar_string_.at(index);
}

inline CigarString::Iterator CigarString::begin() noexcept
{
    return the_cigar_string_.begin();
}

inline CigarString::Iterator CigarString::end() noexcept
{
    return the_cigar_string_.end();
}

inline CigarString::ConstIterator CigarString::begin() const noexcept
{
    return the_cigar_string_.begin();
}

inline CigarString::ConstIterator CigarString::end() const noexcept
{
    return the_cigar_string_.end();
}

inline CigarString::ConstIterator CigarString::cbegin() const noexcept
{
    return the_cigar_string_.cbegin();
}

inline CigarString::ConstIterator CigarString::cend() const noexcept
{
    return the_cigar_string_.cend();
}

inline const CigarString::CigarOperation& CigarString::at(uint_fast32_t n) const
{
    return the_cigar_string_.at(n);
}

inline bool CigarString::empty() const noexcept
{
    return the_cigar_string_.empty();
}

inline uint_fast32_t CigarString::size() const noexcept
{
    return static_cast<uint_fast32_t>(the_cigar_string_.size());
}

inline bool CigarString::is_soft_clipped() const noexcept
{
    return the_cigar_string_[0].get_flag() == 'S';
}

// Defining member and non-member size as the method for calculating
inline uint_fast32_t size(const CigarString& a_cigar_string) noexcept
{
    return a_cigar_string.size();
}

inline std::string to_string(const CigarString::CigarOperation& a_cigar_operation)
{
    return std::to_string(a_cigar_operation.get_size()) + a_cigar_operation.get_flag();
}

inline std::string to_string(const CigarString& a_cigar_string)
{
    std::string result {};
    result.reserve(2 * size(a_cigar_string));
    for (const auto& a_cigar_operation : a_cigar_string) {
        result += to_string(a_cigar_operation);
    }
    return result;
}

inline bool operator==(const CigarString::CigarOperation& lhs, const CigarString::CigarOperation& rhs)
{
    return lhs.get_size() == rhs.get_size() && lhs.get_size() == rhs.get_size();
}

inline bool operator==(const CigarString& lhs, const CigarString& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

inline std::ostream& operator<<(std::ostream& os, const CigarString::CigarOperation& a_cigar_operation)
{
    os << a_cigar_operation.get_size() << a_cigar_operation.get_flag();
    return os;
}

inline std::ostream& operator<<(std::ostream& os, const CigarString& a_cigar_string)
{
    std::copy(std::cbegin(a_cigar_string), std::cend(a_cigar_string),
              std::ostream_iterator<CigarString::CigarOperation>(os));
    return os;
}

#endif
