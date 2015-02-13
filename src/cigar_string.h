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

using std::uint_fast32_t;

class CigarString
{
public:
    using CigarOperation = std::pair<uint_fast32_t, char>;
    using Iterator       = std::vector<CigarOperation>::const_iterator;
    
    CigarString() = delete;
    CigarString(const std::string& a_cigar_string);
    CigarString(std::vector<CigarOperation>&& a_cigar_string);
    
    CigarOperation get_cigar_pair(uint_fast32_t index) const;
    
    Iterator begin() const;
    Iterator end() const;
    
private:
    std::vector<CigarOperation> the_cigar_string_;
};

inline
CigarString::CigarString(const std::string& a_cigar_string) : the_cigar_string_ {}
{
    const uint_fast32_t num_pairs {static_cast<uint_fast32_t>(a_cigar_string.length() / 2)};
    the_cigar_string_.reserve(num_pairs);
    for (uint_fast32_t i {0}; i < num_pairs; ++i) {
        the_cigar_string_.emplace_back(std::make_pair(a_cigar_string[i], a_cigar_string[i + 1]));
    }
}

inline
CigarString::CigarString(std::vector<CigarOperation>&& a_cigar_string)
    : the_cigar_string_ {std::move(a_cigar_string)}
{}

inline
CigarString::CigarOperation CigarString::get_cigar_pair(uint_fast32_t index) const
{
    return the_cigar_string_.at(index);
}

inline
CigarString::Iterator CigarString::begin() const
{
    return std::cbegin(the_cigar_string_);
}

inline
CigarString::Iterator CigarString::end() const
{
    return std::cend(the_cigar_string_);
}

inline
std::string to_string(const CigarString& a_cigar_string)
{
    return std::string {};
}

inline
std::ostream& operator<<(std::ostream& os, const CigarString& a_cigar_string)
{
    os << to_string(a_cigar_string);
    return os;
}

#endif
