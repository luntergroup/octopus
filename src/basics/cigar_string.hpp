//
//  cigar.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_cigar_hpp
#define Octopus_cigar_hpp

#include <string>
#include <cstdint>
#include <iterator>
#include <vector>
#include <numeric>
#include <iosfwd>
#include <functional>

#include <boost/functional/hash.hpp>

#include <interfaces/comparable.hpp>

namespace octopus {

class CigarOperation : public Comparable<CigarOperation> // Comparable so can compare reads
{
public:
    using Size = std::uint_fast32_t;
    using Flag = char;
    
    static constexpr char ALIGNMENT_MATCH {'M'};
    static constexpr char SEQUENCE_MATCH  {'='};
    static constexpr char SUBSTITUTION    {'X'};
    static constexpr char INSERTION       {'I'};
    static constexpr char DELETION        {'D'};
    static constexpr char SOFT_CLIPPED    {'S'};
    static constexpr char HARD_CLIPPED    {'H'};
    static constexpr char PADDING         {'P'};
    static constexpr char SKIPPED         {'N'};
    
    CigarOperation() = default;
    
    explicit CigarOperation(Size size, Flag type) noexcept;
    
    CigarOperation(const CigarOperation&)            = default;
    CigarOperation& operator=(const CigarOperation&) = default;
    CigarOperation(CigarOperation&&)                 = default;
    CigarOperation& operator=(CigarOperation&&)      = default;
    
    ~CigarOperation() = default;
    
    Size size() const noexcept;
    Flag flag() const noexcept;
    
    bool advances_reference() const noexcept;
    bool advances_sequence() const noexcept;
    
private:
    Size size_;
    Flag flag_;
};

bool is_match(const CigarOperation& op) noexcept;
bool is_indel(const CigarOperation& op) noexcept;
bool is_clipping(const CigarOperation& op) noexcept;

// CigarString

using CigarString = std::vector<CigarOperation>;

CigarString parse_cigar(const std::string& cigar);

// Valid if non-empty and all operations are valid
bool is_valid(const CigarString& cigar) noexcept;

// Minimal if all adjacent operations are unique
bool is_minimal(const CigarString& cigar) noexcept;

bool is_front_soft_clipped(const CigarString& cigar) noexcept;

bool is_back_soft_clipped(const CigarString& cigar) noexcept;

bool is_soft_clipped(const CigarString& cigar) noexcept;

std::pair<CigarOperation::Size, CigarOperation::Size> get_soft_clipped_sizes(const CigarString& cigar) noexcept;

template <typename S = CigarOperation::Size>
S clipped_begin(const CigarString& cigar, S unclipped_begin) noexcept
{
    if (is_front_soft_clipped(cigar)) {
        unclipped_begin -= static_cast<S>(cigar.front().size());
    }
    return unclipped_begin;
}

template <typename S = CigarOperation::Size>
S sum_operation_sizes(const CigarString& cigar) noexcept
{
    return std::accumulate(std::cbegin(cigar), std::cend(cigar), S {0},
                           [] (const S curr, const CigarOperation& op) {
                               return curr + op.size();
                           });
}

template <typename S = CigarOperation::Size>
S reference_size(const CigarString& cigar) noexcept
{
    return std::accumulate(std::cbegin(cigar), std::cend(cigar), S {0},
                           [] (const S curr, const CigarOperation& op) {
                               return curr + ((op.advances_reference()) ? op.size() : 0);
                           });
}

template <typename S = CigarOperation::Size>
S sequence_size(const CigarString& cigar) noexcept
{
    return std::accumulate(std::cbegin(cigar), std::cend(cigar), S {0},
                           [] (const S curr, const CigarOperation& op) {
                               return curr + ((op.advances_sequence()) ? op.size() : 0);
                           });
}

template <typename S = CigarOperation::Size>
CigarOperation get_operation_at_sequence_position(const CigarString& cigar, S pos)
{
    auto first = std::cbegin(cigar);
    
    while (pos >= first->size()) {
        ++first;
        pos -= first->size();
    }
    
    return *first;
}

// Relative to both reference and sequence
CigarString splice(const CigarString& cigar, CigarOperation::Size offset,
                   CigarOperation::Size size);

CigarString splice_reference(const CigarString& cigar, CigarOperation::Size offset,
                             CigarOperation::Size size);

CigarString splice_sequence(const CigarString& cigar, CigarOperation::Size offset,
                            CigarOperation::Size size);

inline bool operator==(const CigarOperation& lhs, const CigarOperation& rhs)
{
    return lhs.flag() == rhs.flag() && lhs.size() == rhs.size();
}

inline bool operator<(const CigarOperation& lhs, const CigarOperation& rhs)
{
    return (lhs.flag() == rhs.flag()) ? lhs.size() < rhs.size() :
                                                lhs.flag() < rhs.flag();
}

std::ostream& operator<<(std::ostream& os, const CigarOperation& cigar_operation);

std::ostream& operator<<(std::ostream& os, const CigarString& cigar);

struct CigarHash
{
    std::size_t operator()(const CigarOperation& op) const noexcept;
    std::size_t operator()(const CigarString& cigar) const noexcept;
};

} // namespace octopus

namespace std {
    template <> struct hash<octopus::CigarOperation>
    {
        size_t operator()(const octopus::CigarOperation& op) const noexcept
        {
            return octopus::CigarHash()(op);
        }
    };
    
    template <> struct hash<octopus::CigarString>
    {
        size_t operator()(const octopus::CigarString& cigar) const noexcept
        {
            return octopus::CigarHash()(cigar);
        }
    };
} // namespace std

namespace boost {
    template <>
    struct hash<octopus::CigarOperation> : std::unary_function<octopus::CigarOperation, std::size_t>
    {
        std::size_t operator()(const octopus::CigarOperation& op) const noexcept
        {
            return std::hash<octopus::CigarOperation>()(op);
        }
    };
    
    template <>
    struct hash<octopus::CigarString> : std::unary_function<octopus::CigarString, std::size_t>
    {
        std::size_t operator()(const octopus::CigarString& cigar) const noexcept
        {
            return std::hash<octopus::CigarString>()(cigar);
        }
    };
} // namespace boost

#endif
