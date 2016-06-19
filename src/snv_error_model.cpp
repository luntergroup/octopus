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
    auto extract_repeats(const Haplotype& haplotype)
    {
        return Tandem::find_maximal_repetitions(haplotype.sequence(), 1, 3);
    }
    
    template <typename ForwardIt, typename OutputIt>
    OutputIt count_runs(ForwardIt first, ForwardIt last, OutputIt result,
                        const unsigned max_gap = 1)
    {
        using ValueType = typename std::iterator_traits<ForwardIt>::value_type;
        
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
    
    template <typename C, typename T>
    static auto get_penalty(const C& penalties, const T length)
    {
        return (length <= penalties.size()) ? penalties[length] : penalties.back();
    }
    
    template <typename T, typename C>
    void convert_to_priors(std::vector<T>& run_lengths, const C& penalties)
    {
        std::transform(std::cbegin(run_lengths), std::cend(run_lengths), std::begin(run_lengths),
                       [&penalties] (const auto l) { return get_penalty(penalties, l); } );
    }
    
    constexpr auto base_hash(const char b) noexcept
    {
        using T = std::int8_t;
        switch(b) {
            case 'A': return T {1};
            case 'C': return T {2};
            case 'G': return T {3};
            case 'T': return T {4};
            default: return T {5};
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
    } // namespace
    
    void SnvErrorModel::evaluate(const Haplotype& haplotype,
                                     PenaltyVector& forward_snv_priors,
                                     PenaltyVector& reverse_snv_priors) const
    {
        const auto repeats = extract_repeats(haplotype);
        
        std::vector<std::int8_t> repeat_mask(sequence_size(haplotype), 0);
        
        for (const auto& repeat : repeats) {
            std::fill_n(std::next(std::begin(repeat_mask), repeat.pos), repeat.length,
                        repeat_hash(haplotype, repeat));
        }
        
        forward_snv_priors.resize(sequence_size(haplotype));
        reverse_snv_priors.resize(sequence_size(haplotype));
        
        count_runs(std::cbegin(repeat_mask), std::cend(repeat_mask), std::begin(forward_snv_priors));
        count_runs(std::crbegin(repeat_mask), std::crend(repeat_mask), std::rbegin(reverse_snv_priors));
        
//        ::debug::print_variant_alleles(haplotype); std::cout << '\n';
//        std::cout << haplotype.sequence() << std::endl;
//        std::copy(std::cbegin(repeat_mask), std::cend(repeat_mask),
//                  std::ostream_iterator<int>(std::cout));
//        std::cout << std::endl;
//        std::transform(std::cbegin(forward_snv_priors), std::cend(forward_snv_priors),
//                       std::ostreambuf_iterator<char>(std::cout),
//                       [] (const auto i) -> char { return std::min(i + 48, 126); });
//        std::cout << std::endl;
//        std::transform(std::cbegin(reverse_snv_priors), std::cend(reverse_snv_priors),
//                       std::ostreambuf_iterator<char>(std::cout),
//                       [] (const auto i) -> char { return std::min(i + 48, 126); });
//        std::cout << std::endl;
//        exit(0);
        
        convert_to_priors(forward_snv_priors, Max_qualities_);
        convert_to_priors(reverse_snv_priors, Max_qualities_);
    }
} // namespace Octopus
