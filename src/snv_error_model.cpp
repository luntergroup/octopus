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

#include "haplotype.hpp"
#include "tandem.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
    constexpr decltype(SnvErrorModel::Homopolymer_errors_) SnvErrorModel::Homopolymer_errors_;
    
    namespace
    {
        auto extract_repeats(const Haplotype& haplotype)
        {
            return Tandem::find_maximal_repetitions(haplotype.sequence(), 1, 3);
        }
    }
    
    template <typename C, typename T>
    static auto get_penalty(const C& penalties, const T length)
    {
        return (length <= penalties.size()) ? penalties[length] : penalties.back();
    }
    
    template <typename ForwardIt, typename OutputIt>
    OutputIt count_runs(ForwardIt first, ForwardIt last, OutputIt result,
                        const unsigned max_gap = 1)
    {
        using ValueType = typename std::iterator_traits<ForwardIt>::value_type;
        
        ValueType count {0};
        unsigned gap {0};
        
        return std::transform(first, last, result,
                              [&count, &gap, max_gap] (const auto x) -> ValueType {
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
                                  } else {
                                      gap = 0;
                                      count += x;
                                  }
                                  return 0;
                              });
    }
    
    template <typename T, typename C>
    void convert_to_penalties(std::vector<T>& run_lengths, const C& penalties)
    {
        std::transform(std::cbegin(run_lengths), std::cend(run_lengths), std::begin(run_lengths),
                       [&penalties] (const auto l) { return get_penalty(penalties, l); } );
    }
    
    void SnvErrorModel::evaluate(const Haplotype& haplotype,
                                     PenaltyVector& forward_snv_priors,
                                     PenaltyVector& reverse_snv_priors) const
    {
        const auto repeats = extract_repeats(haplotype);
        
        std::vector<std::int8_t> repeat_mask(sequence_size(haplotype), 0);
        
        for (const auto& repeat : repeats) {
            std::fill_n(std::next(std::begin(repeat_mask), repeat.pos), repeat.length, 1);
        }
        
        forward_snv_priors.resize(sequence_size(haplotype));
        reverse_snv_priors.resize(sequence_size(haplotype));
        
        count_runs(std::cbegin(repeat_mask), std::cend(repeat_mask), std::begin(forward_snv_priors));
        count_runs(std::crbegin(repeat_mask), std::crend(repeat_mask), std::rbegin(reverse_snv_priors));
        
        ::debug::print_variant_alleles(haplotype); std::cout << '\n';
        std::cout << haplotype.sequence() << std::endl;
        std::copy(std::cbegin(repeat_mask), std::cend(repeat_mask), std::ostream_iterator<int>(std::cout));
        std::cout << std::endl;
        std::transform(std::cbegin(forward_snv_priors), std::cend(forward_snv_priors),
                       std::ostreambuf_iterator<char>(std::cout),
                       [] (const auto i) -> char { return std::min(i + 48, 126); });
        std::cout << std::endl;
        std::transform(std::cbegin(reverse_snv_priors), std::cend(reverse_snv_priors),
                       std::ostreambuf_iterator<char>(std::cout),
                       [] (const auto i) -> char { return std::min(i + 48, 126); });
        std::cout << std::endl;
        //exit(0);
        
        convert_to_penalties(forward_snv_priors, Homopolymer_errors_);
        convert_to_penalties(reverse_snv_priors, Homopolymer_errors_);
    }
} // namespace Octopus
