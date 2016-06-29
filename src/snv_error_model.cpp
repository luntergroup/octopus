//
//  snv_error_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/06/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "snv_error_model.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>

#include "haplotype.hpp"
#include "tandem.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
    constexpr decltype(SnvErrorModel::Max_qualities_) SnvErrorModel::Max_qualities_;
    
    namespace
    {
    auto extract_repeats(const Haplotype& haplotype, const unsigned max_period)
    {
        return Tandem::find_maximal_repetitions(haplotype.sequence(), 1, max_period);
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
                              [&count, &gap, &prev, max_gap] (const auto x) -> ValueType {
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
        switch(b) {
            case 'A': return T {1};
            case 'C': return T {2};
            case 'G': return T {3};
            case 'T': return T {4};
            default:  return T {5};
        }
    }
    
    auto repeat_hash(const Haplotype& haplotype, const Tandem::StringRun& repeat)
    {
        const auto& sequence = haplotype.sequence();
        const auto first = std::next(std::begin(sequence), repeat.pos);
        const auto last  = std::next(first, repeat.period);
        return std::accumulate(first, last, std::int8_t {0},
                               [] (const auto& curr, const auto b) {
                                   return curr + base_hash(b);
                               });
    }
    
    template <typename C, typename T>
    static auto get_penalty(const C& penalties, const T length)
    {
        return (length <= penalties.size()) ? penalties[length] : penalties.back();
    }
    
    template <typename T1, typename T2, typename C>
    void set_priors(const std::vector<T1>& run_lengths, std::vector<T2>& result, const C& penalties)
    {
        std::transform(std::cbegin(run_lengths), std::cend(run_lengths), std::cbegin(result),
                       std::begin(result),
                       [&penalties] (const auto l, const auto curr) {
                           return std::min(get_penalty(penalties, l), curr);
                       } );
    }
    } // namespace
    
    void SnvErrorModel::evaluate(const Haplotype& haplotype,
                                 PenaltyVector& forward_snv_priors,
                                 PenaltyVector& reverse_snv_priors) const
    {
        constexpr auto Max_period = Max_qualities_.size();
        
        const auto repeats = extract_repeats(haplotype, Max_period);
        
        const auto num_bases = sequence_size(haplotype);
        
        std::array<std::vector<std::int8_t>, Max_period> repeat_masks {};
        
        repeat_masks.fill(std::vector<std::int8_t>(num_bases, 0));
        
        for (const auto& repeat : repeats) {
            std::fill_n(std::next(std::begin(repeat_masks[repeat.period - 1]), repeat.pos), repeat.length,
                        repeat_hash(haplotype, repeat));
        }
        
        forward_snv_priors.assign(num_bases, Max_qualities_.front().front());
        reverse_snv_priors.assign(num_bases, Max_qualities_.front().front());
        
        std::vector<unsigned> runs(num_bases);
        
        for (std::size_t i {0}; i < Max_period; ++i) {
            const auto& repeat_mask = repeat_masks[i];
            count_runs(std::cbegin(repeat_mask), std::cend(repeat_mask), std::begin(runs));
            set_priors(runs, forward_snv_priors, Max_qualities_[i]);
            count_runs(std::crbegin(repeat_mask), std::crend(repeat_mask), std::rbegin(runs));
            set_priors(runs, reverse_snv_priors, Max_qualities_[i]);
        }
    }
} // namespace Octopus
