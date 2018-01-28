// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "x10_snv_error_model.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>

#include <tandem/tandem.hpp>

#include <core/types/haplotype.hpp>

namespace octopus {

constexpr decltype(X10SnvErrorModel::maxQualities_) X10SnvErrorModel::maxQualities_;

std::unique_ptr<SnvErrorModel> X10SnvErrorModel::do_clone() const
{
    return std::make_unique<X10SnvErrorModel>(*this);
}

namespace {

auto extract_repeats(const Haplotype& haplotype, const unsigned max_period)
{
    return tandem::extract_exact_tandem_repeats(haplotype.sequence(), 1, max_period);
}

template <typename ForwardIt, typename OutputIt>
OutputIt count_runs(ForwardIt first, ForwardIt last, OutputIt result,
                    const unsigned max_gap = 4)
{
    using ValueType = typename std::iterator_traits<OutputIt>::value_type;
    if (first == last) return result;
    auto prev = *first;
    auto count = (prev > 0) ? 1 : 0;
    unsigned gap {0};
    *result++ = 0;
    return std::transform(std::next(first), last, result,
                          [&count, &gap, &prev, max_gap](const auto x) -> ValueType {
                              if (x == 0) {
                                  ++gap;
                                  if (count > 0) {
                                      if (gap == 1) {
                                          if (max_gap >= 1) {
                                              return count;
                                          } else {
                                              const auto tmp = count;
                                              count = 0;
                                              return tmp;
                                          }
                                      } else if (gap > max_gap) {
                                          count = 0;
                                      }
                                  }
                              } else if (prev == x) {
                                  gap = 0;
                                  ++count;
                              } else {
                                  prev = x;
                                  const auto tmp = count;
                                  count = 1;
                                  return tmp;
                              }
                              return 0;
                          });
}

constexpr auto base_hash(const char b) noexcept
{
    using T = std::int8_t;
    switch (b) {
    case 'A':
        return T {1};
    case 'C':
        return T {2};
    case 'G':
        return T {3};
    case 'T':
        return T {4};
    default:
        return T {5};
    }
}

auto repeat_hash(const Haplotype& haplotype, const tandem::Repeat& repeat) noexcept
{
    const auto& sequence = haplotype.sequence();
    const auto first = std::next(std::begin(sequence), repeat.pos);
    const auto last = std::next(first, repeat.period);
    return std::accumulate(first, last, std::int8_t {0}, [] (const auto& curr, const auto b) { return curr + base_hash(b); });
}

template <typename C, typename T>
static auto get_penalty(const C& penalties, const T length)
{
    return length < penalties.size() ? penalties[length] : penalties.back();
}

template <typename T1, typename T2, typename C>
void set_priors(const std::vector<T1>& run_lengths, std::vector<T2>& result, const C& penalties)
{
    std::transform(std::cbegin(run_lengths), std::cend(run_lengths), std::cbegin(result),
                   std::begin(result),
                   [&penalties](const auto l, const auto curr) {
                       return std::min(get_penalty(penalties, l), curr);
                   });
}

auto make_substitution_mask(const Haplotype& haplotype)
{
    const auto cigar = haplotype.cigar();
    std::vector<bool> result(sequence_size(haplotype));
    auto mask_itr = std::begin(result);
    for (const auto& op : cigar) {
        if (advances_sequence(op)) {
            mask_itr = std::fill_n(mask_itr, op.size(), op.flag() == CigarOperation::Flag::substitution);
        }
    }
    return result;
}

} // namespace

void X10SnvErrorModel::do_evaluate(const Haplotype& haplotype,
                                   MutationVector& forward_snv_mask, PenaltyVector& forward_snv_priors,
                                   MutationVector& reverse_snv_mask, PenaltyVector& reverse_snv_priors) const
{
    using std::cbegin; using std::cend; using std::crbegin; using std::crend;
    using std::begin; using std::end; using std::rbegin; using std::next;
    constexpr auto Max_period = maxQualities_.size();
    const auto repeats = extract_repeats(haplotype, Max_period);
    const auto num_bases = sequence_size(haplotype);
    std::array<std::vector<std::int8_t>, Max_period> repeat_masks {};
    repeat_masks.fill(std::vector<std::int8_t>(num_bases, 0));
    for (const auto& repeat : repeats) {
        std::fill_n(next(begin(repeat_masks[repeat.period - 1]), repeat.pos), repeat.length, repeat_hash(haplotype, repeat));
    }
    const auto max_quality = maxQualities_.front().front();
    forward_snv_priors.assign(num_bases, max_quality);
    reverse_snv_priors.assign(num_bases, max_quality);
    std::vector<unsigned> runs(num_bases);
    constexpr unsigned Max_homopolymer_gap {3};
    for (std::int8_t i {1}; i < 5; ++i) {
        auto homopolymers = repeat_masks[0];
        std::replace_if(begin(homopolymers), end(homopolymers), [=] (auto x) { return x != i; }, 0);
        count_runs(cbegin(homopolymers), cend(homopolymers), begin(runs), Max_homopolymer_gap);
        set_priors(runs, forward_snv_priors, maxQualities_.front());
        count_runs(crbegin(homopolymers), crend(homopolymers), rbegin(runs), Max_homopolymer_gap);
        set_priors(runs, reverse_snv_priors, maxQualities_.front());
    }
    for (std::size_t i {1}; i < Max_period; ++i) {
        const auto max_gap = Max_homopolymer_gap + i;
        const auto& repeat_mask = repeat_masks[i];
        count_runs(cbegin(repeat_mask), cend(repeat_mask), begin(runs), max_gap);
        set_priors(runs, forward_snv_priors, maxQualities_[i]);
        count_runs(crbegin(repeat_mask), crend(repeat_mask), rbegin(runs), max_gap);
        set_priors(runs, reverse_snv_priors, maxQualities_[i]);
    }
    const auto substitution_mask = make_substitution_mask(haplotype);
    std::transform(std::cbegin(forward_snv_priors), std::cend(forward_snv_priors), std::cbegin(substitution_mask),
                   std::begin(forward_snv_priors), [=] (auto q, auto b) { return !b ? q : max_quality; });
    std::transform(std::cbegin(reverse_snv_priors), std::cend(reverse_snv_priors), std::cbegin(substitution_mask),
                   std::begin(reverse_snv_priors), [=] (auto q, auto b) { return !b ? q : max_quality; });
    const auto& sequence = haplotype.sequence();
    forward_snv_mask.resize(num_bases);
    std::rotate_copy(crbegin(sequence), next(crbegin(sequence)), crend(sequence), rbegin(forward_snv_mask));
    reverse_snv_mask.resize(num_bases);
    std::rotate_copy(cbegin(sequence), next(cbegin(sequence)), cend(sequence), begin(reverse_snv_mask));
}

} // namespace octopus
