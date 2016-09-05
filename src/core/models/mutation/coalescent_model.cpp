// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coalescent_model.hpp"

#include <stdexcept>

#include "tandem/tandem.hpp"

namespace octopus {

auto percent_of_bases_in_repeat(const Haplotype& haplotype)
{
    const auto repeats = tandem::find_maximal_repetitions(haplotype.sequence(), 1, 6);
    
    if (repeats.empty()) return 0.0;
    
    std::vector<unsigned> repeat_counts(sequence_size(haplotype), 0);
    
    for (const auto& repeat : repeats) {
        const auto it1 = std::next(std::begin(repeat_counts), repeat.pos);
        const auto it2 = std::next(it1, repeat.length);
        std::transform(it1, it2, it1, [] (const auto c) { return c + 1; });
    }
    
    const auto c = std::count_if(std::cbegin(repeat_counts), std::cend(repeat_counts),
                                 [] (const auto c) { return c > 0; });
    return static_cast<double>(c) / repeat_counts.size();
}

auto calculate_base_indel_heterozygosities(const Haplotype& haplotype,
                                           const double base_indel_heterozygosity)
{
    std::vector<double> result(sequence_size(haplotype), base_indel_heterozygosity);
    
    const auto repeats = tandem::find_maximal_repetitions(haplotype.sequence(), 1, 3);
    
    for (const auto& repeat : repeats) {
        const auto it1 = std::next(std::begin(result), repeat.pos);
        const auto it2 = std::next(it1, repeat.length);
        const auto n = repeat.length / repeat.period;
        // TODO: implement a proper model for this
        const auto t = std::min(base_indel_heterozygosity * std::pow(n, 2.6), 1.0);
        std::transform(it1, it2, it1, [t] (const auto h) { return std::max(h, t); });
    }
    
    return result;
}

CoalescentModel::CoalescentModel(Haplotype reference,
                                 Parameters params,
                                 unsigned max_haplotypes)
: reference_ {std::move(reference)}
, reference_base_indel_heterozygosities_ {}
, params_ {params}
{
    if (params_.snp_heterozygosity <= 0 || params_.indel_heterozygosity <= 0) {
        throw std::domain_error {"CoalescentModel: snp and indel heterozygosity must be > 0"};
    }
    
    reference_base_indel_heterozygosities_ = calculate_base_indel_heterozygosities(reference_, params_.indel_heterozygosity);
    
    site_buffer1_.reserve(128);
    site_buffer2_.reserve(128);
    difference_cache_.reserve(max_haplotypes);
    difference_cache_.emplace(std::piecewise_construct,
                              std::forward_as_tuple(reference_),
                              std::forward_as_tuple());
    result_cache_.reserve(max_haplotypes);
}

void CoalescentModel::set_reference(Haplotype reference)
{
    reference_ = std::move(reference);
    difference_cache_.clear();
    difference_cache_.emplace(std::piecewise_construct,
                              std::forward_as_tuple(reference_),
                              std::forward_as_tuple());
}

} // namespace octopus
