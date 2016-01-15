//
//  cigar_string.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_cigar_string_hpp
#define Octopus_cigar_string_hpp

#include <string>
#include <cstdint>
#include <iterator>
#include <vector>
#include <ostream>
#include <numeric>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "comparable.hpp"

class CigarOperation : public Comparable<CigarOperation> // Comparable so can compare reads
{
public:
    using SizeType = std::uint_fast32_t;
    
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
    explicit CigarOperation(SizeType size, char type) noexcept;
    ~CigarOperation() = default;
    
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

boost::optional<CigarString> parse_cigar_string(const std::string& cigar_string);

bool is_valid_cigar(const CigarString& cigar) noexcept;

bool is_minimal_cigar(const CigarString& cigar) noexcept;

bool is_front_soft_clipped(const CigarString& cigar_string) noexcept;

bool is_back_soft_clipped(const CigarString& cigar_string) noexcept;

bool is_soft_clipped(const CigarString& cigar_string) noexcept;

std::pair<CigarOperation::SizeType, CigarOperation::SizeType>
get_soft_clipped_sizes(const CigarString& cigar_string) noexcept;

template <typename SizeType = CigarOperation::SizeType>
SizeType soft_clipped_read_begin(const CigarString& cigar_string, SizeType hard_clipped_begin) noexcept
{
    if (is_front_soft_clipped(cigar_string)) {
        hard_clipped_begin -= static_cast<SizeType>(cigar_string.front().get_size());
    }
    return hard_clipped_begin;
}

template <typename SizeType = CigarOperation::SizeType>
SizeType operations_size(const CigarString& cigar_string) noexcept
{
    return std::accumulate(std::cbegin(cigar_string), std::cend(cigar_string), SizeType {0},
                           [] (const SizeType curr, const CigarOperation& op) {
                               return curr + op.get_size();
                           });
}

template <typename SizeType = CigarOperation::SizeType>
SizeType reference_size(const CigarString& cigar_string) noexcept
{
    return std::accumulate(std::cbegin(cigar_string), std::cend(cigar_string), SizeType {0},
                           [] (const SizeType curr, const CigarOperation& op) {
                               return curr + ((op.advances_reference()) ? op.get_size() : 0);
                           });
}

template <typename SizeType = CigarOperation::SizeType>
SizeType sequence_size(const CigarString& cigar_string) noexcept
{
    return std::accumulate(std::cbegin(cigar_string), std::cend(cigar_string), SizeType {0},
                           [] (const SizeType curr, const CigarOperation& op) {
                               return curr + ((op.advances_sequence()) ? op.get_size() : 0);
                           });
}

template <typename SizeType = CigarOperation::SizeType>
CigarOperation get_operation_at_sequence_position(const CigarString& cigar_string,
                                                  const SizeType sequence_pos)
{
    auto first = std::cbegin(cigar_string);
    
    while (sequence_pos >= first->get_size()) {
        ++first;
        sequence_pos -= first->get_size();
    }
    
    return *first;
}

CigarString splice(const CigarString& cigar_string, CigarOperation::SizeType offset,
                   CigarOperation::SizeType size);

CigarString splice(const CigarString& cigar_string, CigarOperation::SizeType size);

CigarString splice_reference(const CigarString& cigar_string, CigarOperation::SizeType offset,
                             CigarOperation::SizeType size);

CigarString splice_reference(const CigarString& cigar_string, CigarOperation::SizeType size);

CigarString splice_sequence(const CigarString& cigar_string, CigarOperation::SizeType offset,
                            CigarOperation::SizeType size);

CigarString splice_sequence(const CigarString& cigar_string, CigarOperation::SizeType size);

inline bool operator==(const CigarOperation& lhs, const CigarOperation& rhs)
{
    return lhs.get_flag() == rhs.get_flag() && lhs.get_size() == rhs.get_size();
}

inline bool operator<(const CigarOperation& lhs, const CigarOperation& rhs)
{
    return (lhs.get_flag() == rhs.get_flag()) ? lhs.get_size() < rhs.get_size() :
                                                lhs.get_flag() < rhs.get_flag();
}

std::ostream& operator<<(std::ostream& os, const CigarOperation& cigar_operation);

std::ostream& operator<<(std::ostream& os, const CigarString& cigar_string);

namespace std {
    template <> struct hash<CigarOperation>
    {
        size_t operator()(const CigarOperation& op) const
        {
            using boost::hash_combine;
            size_t result {};
            hash_combine(result, op.get_flag());
            hash_combine(result, op.get_size());
            return result;
        }
    };
    
    template <> struct hash<CigarString>
    {
        size_t operator()(const CigarString& cigar) const
        {
            return boost::hash_range(std::cbegin(cigar), std::cend(cigar));
        }
    };
} // namespace std

namespace boost {
    template <> struct hash<CigarOperation> : std::unary_function<CigarOperation, std::size_t>
    {
        std::size_t operator()(const CigarOperation& op) const
        {
            return std::hash<CigarOperation>()(op);
        }
    };
} // namespace boost

#endif
