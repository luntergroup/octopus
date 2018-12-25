// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef sequence_utils_hpp
#define sequence_utils_hpp

#include <vector>
#include <array>
#include <unordered_map>
#include <map>
#include <cstddef>
#include <iterator>
#include <algorithm>
#include <functional>
#include <random>

#include <boost/optional.hpp>

#include "tandem/tandem.hpp"

#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"

namespace octopus { namespace utils {

namespace detail {

static constexpr std::array<char, 4> dnaBases {'A', 'C', 'G', 'T'};
static constexpr std::array<char, 4> rnaBases {'A', 'C', 'G', 'U'};

static const std::array<std::vector<char>, 128> iupac_symbols
{{
 {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},
 {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},
 {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},
 {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {},
 {},
 {'A'}, // Adenine
 {'C', 'G', 'T', 'U'}, // B (not A)
 {'C'}, // Cytosine
 {'A', 'G', 'T', 'U'}, // D (not C)
 {}, {},
 {'G'}, // Guanine
 {'A', 'C', 'T', 'U'}, // H (not G)
 {}, {},
 {'G', 'T', 'U'}, // Ketones
 {},
 {'A', 'C'}, // aMino groups
 {'A', 'C', 'G', 'T', 'U'}, // Nucleic acid
 {}, {}, {},
 {'A', 'G'}, // puRine
 {'C', 'G'}, // Strong interaction
 {'T'}, // Thymine
 {'U'}, // Uracil
 {'A', 'C', 'G'}, // V (not T/U)
 {'A', 'T', 'U'}, // Weak interaction
 {},
 {'C', 'T', 'U'}, // pYrimidines
 {}, {}, {}, {}, {}, {}, {},
 {'a'}, // adenine
 {'c', 'g', 't', 'u'}, // B (not a)
 {'c'}, // cytosine
 {'a', 'g', 't', 'u'}, // D (not c)
 {}, {},
 {'g'}, // guanine
 {'a', 'c', 't', 'u'}, // H (not g)
 {}, {},
 {'g', 't', 'u'}, // Ketones
 {},
 {'a', 'c'}, // aMino groups
 {'a', 'c', 'g', 't', 'u'}, // Nucleic acid
 {}, {}, {},
 {'a', 'g'}, // puRine
 {'c', 'g'}, // Strong interaction
 {'t'}, // thymine
 {'u'}, // uracil
 {'a', 'c', 'g'}, // v (not t/u)
 {'a', 't', 'u'}, // Weak interaction
 {},
 {'c', 't', 'u'}, // pYrimidines
 {}, {}, {}, {}, {}, {}
 }};

inline constexpr char capitalise_base(const char base) noexcept
{
    switch (base) {
        case 'a': return 'A';
        case 'b': return 'B';
        case 'c': return 'C';
        case 'd': return 'C';
        case 'g': return 'G';
        case 'h': return 'H';
        case 'k': return 'K';
        case 't': return 'T';
        case 'u': return 'U';
        case 'm': return 'M';
        case 'n': return 'N';
        case 'r': return 'R';
        case 's': return 'S';
        case 'v': return 'V';
        case 'w': return 'W';
        case 'y': return 'Y';
        default : return base;
    }
}

inline bool is_dna_nucleotide(const char b) noexcept
{
    return b == 'A' || b == 'C' || b == 'G' || b == 'T';
}

inline bool is_rna_nucleotide(const char b) noexcept
{
    return b == 'A' || b == 'C' || b == 'G' || b == 'U';
}

inline bool is_dna_or_rna_nucleotide(const char b) noexcept
{
    return b == 'A' || b == 'C' || b == 'G' || b == 'T' || b == 'U';
}

inline bool is_iupac_ambiguous_base(const char base) noexcept
{
    return iupac_symbols[base].size() > 1;
}

inline char disambiguate_iupac_base(const char base)
{
    return iupac_symbols[base].front();
}

} // namespace detail

template <typename SequenceType>
bool has_ns(const SequenceType& sequence) noexcept
{
    return std::find(std::cbegin(sequence), std::cend(sequence), 'N') != std::cend(sequence);
}

template <typename SequenceType>
bool is_dna(const SequenceType& sequence) noexcept
{
    return std::all_of(std::cbegin(sequence), std::cend(sequence),
                       [] (const auto c) noexcept { return detail::is_dna_nucleotide(c) || c == 'N'; });
}

template <typename SequenceType>
bool is_canonical_dna(const SequenceType& sequence) noexcept
{
    return std::all_of(std::cbegin(sequence), std::cend(sequence), detail::is_dna_nucleotide);
}

template <typename SequenceType>
bool is_rna(const SequenceType& sequence) noexcept
{
    return std::all_of(std::cbegin(sequence), std::cend(sequence),
                       [] (const auto c) noexcept { return detail::is_rna_nucleotide(c) || c == 'N'; });
}

template <typename SequenceType>
bool is_canonical_rna(const SequenceType& sequence) noexcept
{
    return std::all_of(std::cbegin(sequence), std::cend(sequence), detail::is_rna_nucleotide);
}

template <typename SequenceType>
bool is_dna_rna_ambiguous(const SequenceType& sequence) // i.e. is_dna(sequence) && is_rna(sequence) == true
{
    return sequence.find_first_of("TU") != SequenceType::npos;
}

template <typename SequenceType>
bool is_canonical_dna_or_rna(const SequenceType& sequence) noexcept
{
    return std::all_of(std::cbegin(sequence), std::cend(sequence), detail::is_dna_or_rna_nucleotide);
}

template <typename SequenceType>
void transcribe(SequenceType& dna_sequence)
{
    std::replace(std::begin(dna_sequence), std::end(dna_sequence), 'T', 'U');
}

template <typename SequenceType>
SequenceType transcribe_copy(SequenceType dna_sequence)
{
    std::replace(std::begin(dna_sequence), std::end(dna_sequence), 'T', 'U');
    return dna_sequence;
}

template <typename SequenceType>
void reverse_transcribe(SequenceType& rna_sequence)
{
    std::replace(std::begin(rna_sequence), std::end(rna_sequence), 'U', 'T');
}

template <typename SequenceType>
SequenceType reverse_transcribe_copy(SequenceType rna_sequence)
{
    std::replace(std::begin(rna_sequence), std::end(rna_sequence), 'U', 'T');
    return rna_sequence;
}

template <typename SequenceType>
void capitalise(SequenceType& sequence)
{
    std::transform(std::begin(sequence), std::end(sequence), std::begin(sequence),
                   [] (auto base) { return detail::capitalise_base(base); });
}

template <typename SequenceType>
SequenceType capitalise_copy(const SequenceType& sequence)
{
    auto result = sequence;
    capitalise(result);
    return result;
}

namespace detail {

template <typename Container>
typename Container::value_type random_member(const Container& values)
{
    static std::default_random_engine generator {};
    if (values.empty()) throw std::runtime_error {"cannot sample from empty container"};
    if (values.size() == 1) return *std::cbegin(values);
    std::uniform_int_distribution<size_t> distribution {0, values.size() - 1};
    return *std::next(std::cbegin(values), distribution(generator));
}

} // namespace detail

template <typename SequenceType>
void disambiguate_iupac_bases(SequenceType& sequence)
{
    std::transform(std::begin(sequence), std::end(sequence), std::begin(sequence),
                   [] (auto base) { return detail::disambiguate_iupac_base(base); });
}

template <typename SequenceType>
SequenceType disambiguate_iupac_bases_copy(const SequenceType& sequence)
{
    auto result = sequence;
    disambiguate_iupac_bases(result);
    return result;
}

namespace detail {

static constexpr std::array<char, 128> complementTable
{
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
    4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
    4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 78, 4,
    4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4,  4,
    4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4,  4,
    4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4,  4
};

} // namespace detail

static inline constexpr char complement(const char base) noexcept
{
    return detail::complementTable[base];
}

template <typename BidirIt, typename OutputIt>
OutputIt reverse_complement_copy(BidirIt first, BidirIt last, OutputIt result)
{
    return std::transform(std::make_reverse_iterator(last), std::make_reverse_iterator(first),
                          result, [] (const char base) { return complement(base); });
}

template <typename SequenceType>
SequenceType reverse_complement_copy(const SequenceType& sequence)
{
    SequenceType result {};
    result.resize(sequence.size());
    reverse_complement_copy(std::cbegin(sequence), std::cend(sequence), std::begin(result));
    return result;
}

namespace detail {

static inline void complement_swap(char& lhs, char& rhs) noexcept
{
    const auto tmp = complement(lhs);
    lhs = complement(rhs);
    rhs = tmp;
}

template <typename BidirIt>
void reverse_complement(BidirIt first, BidirIt last, std::bidirectional_iterator_tag)
{
    while (first != last)
    {
        if (first == --last) {
            *first = complement(*first);
            break;
        }
        complement_swap(*first, *last);
        ++first;
    }
}

template <typename BidirIt>
void reverse_complement(BidirIt first, BidirIt last, std::random_access_iterator_tag)
{
    if (first != last)
        for (; first < --last; ++first)
            complement_swap(*first, *last);
    if (first == last)
        *first = complement(*first);
}

} // namespace detail

template <typename BidirIt>
void reverse_complement(BidirIt first, BidirIt last)
{
    detail::reverse_complement(first, last, typename std::iterator_traits<BidirIt>::iterator_category {});
}

template <typename SequenceType>
void reverse_complement(SequenceType& sequence)
{
    reverse_complement(std::begin(sequence), std::end(sequence));
}

template <typename InputIt, typename BidirIt>
bool are_reverse_complements(InputIt first1, InputIt last1, BidirIt first2, BidirIt last2)
{
    return std::equal(first1, last1, std::make_reverse_iterator(last2), std::make_reverse_iterator(first2),
                      [] (const char lhs, const char rhs) noexcept { return lhs == complement(rhs); });
}

template <typename SeqType1, typename SeqType2>
bool are_reverse_complements(const SeqType1& lhs, const SeqType2& rhs)
{
    return are_reverse_complements(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::cend(rhs));
}

namespace detail {

template <typename RandomIt>
bool is_palindromic(RandomIt first, RandomIt last, std::random_access_iterator_tag)
{
    if (first == last) return true;
    
    const auto size = std::distance(first, last);
    
    if (size % 2 != 0) return false;
    
    return std::equal(first, std::next(first, size / 2), std::make_reverse_iterator(last),
                      [] (const char lhs, const char rhs) {
                          return lhs == complement(rhs);
                      });
}

template <typename BidirIt>
bool is_palindromic(BidirIt first, BidirIt last, std::bidirectional_iterator_tag)
{
    if (first == last)   return true;
    if (first == --last) return false;
    
    for (; first != last; ++first, --last) {
        if (*first != complement(*last)) return false;
        if (std::next(first) == last) return true;
    }
    
    return false;
}

} // namespace detail

template <typename BidirIt>
bool is_palindromic(BidirIt first, BidirIt last)
{
    return detail::is_palindromic(first, last, typename std::iterator_traits<BidirIt>::iterator_category {});
}

template <typename SequenceType>
bool is_palindromic(const SequenceType& sequence)
{
    return is_palindromic(std::cbegin(sequence), std::cend(sequence));
}

template <typename SequenceType>
std::unordered_map<char, size_t> count_bases(const SequenceType& sequence)
{
    std::unordered_map<char, size_t> result {};
    result.reserve(5); // 4 bases + N
    
    for (char base : sequence) {
        ++result[base];
    }
    
    return result;
}

template <typename SequenceType>
double gc_content(const SequenceType& sequence) noexcept
{
    const auto gc_count = std::count_if(std::cbegin(sequence), std::cend(sequence),
                                        [] (const char base) { return base == 'G' || base == 'C'; });
    return static_cast<double>(gc_count) / sequence.size();
}

template <typename ForwardIt>
bool is_homopolymer(ForwardIt first, ForwardIt last) noexcept
{
    return std::distance(first, last) > 0 && detail::is_dna_or_rna_nucleotide(*first)
           && std::adjacent_find(first, last, std::not_equal_to<> {}) == last;
}

template <typename SequenceType>
bool is_homopolymer(const SequenceType& sequence) noexcept
{
    return is_homopolymer(std::cbegin(sequence), std::cend(sequence));
}

template <typename Sequence>
bool is_tandem_repeat(const Sequence& sequence, const unsigned period) noexcept
{
    if (period == 0 || sequence.empty() || period > sequence.size() / 2) {
        return false;
    }
    if (period == 1) {
        return is_homopolymer(sequence);
    }
    const auto num_full_periods = sequence.size() / period;
    const auto first_unit_begin_itr = std::cbegin(sequence);
    const auto first_unit_end_itr   = std::next(first_unit_begin_itr, period);
    for (unsigned i {1}; i < num_full_periods; ++i) {
        if (!std::equal(first_unit_begin_itr, first_unit_end_itr, std::next(first_unit_begin_itr, i * period))) {
            return false;
        }
    }
    const auto last_unit_begin_itr = std::next(first_unit_begin_itr, num_full_periods * period);
    if (last_unit_begin_itr != std::cend(sequence)) {
        return std::equal(last_unit_begin_itr, std::cend(sequence), first_unit_begin_itr);
    } else {
        return true;
    }
}

template <typename Sequence>
boost::optional<unsigned> find_tandem_repeat_period(const Sequence& sequence) noexcept
{
    for (unsigned period {1}; period < sequence.size() / 2; ++period) {
        if (is_tandem_repeat(sequence, period)) return period;
    }
    return boost::none;
}

template <typename Sequence>
bool is_tandem_repeat(const Sequence& sequence) noexcept
{
    return find_tandem_repeat_period(sequence);
}

} // namespace utils
} // namespace octopus

#endif
