// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "repeat_based_indel_error_model.hpp"

#include <algorithm>
#include <iterator>

#include "tandem/tandem.hpp"

namespace octopus {

namespace {

auto extract_repeats(const Haplotype& haplotype)
{
    return tandem::extract_exact_tandem_repeats(haplotype.sequence(), 1, 3);
}

void set_motif(const Haplotype& haplotype, const tandem::Repeat& repeat, Haplotype::NucleotideSequence& result)
{
    const auto motif_itr = std::next(std::cbegin(haplotype.sequence()), repeat.pos);
    result.assign(motif_itr, std::next(motif_itr, repeat.period));
}

template <typename FordwardIt, typename Tp>
auto fill_if_less(FordwardIt first, FordwardIt last, const Tp& value)
{
    return std::transform(first, last, first, [&] (const auto& x) { return std::min(x, value); });
}

template <typename FordwardIt, typename Tp>
auto fill_n_if_less(FordwardIt first, std::size_t n, const Tp& value)
{
    return fill_if_less(first, std::next(first, n), value);
}

} // namespace

void RepeatBasedIndelErrorModel::do_set_penalties(const Haplotype& haplotype, PenaltyVector& gap_open_penalities, PenaltyType& gap_extend_penalty) const
{
    gap_open_penalities.assign(sequence_size(haplotype), get_default_open_penalty());
    const auto repeats = extract_repeats(haplotype);
    if (!repeats.empty()) {
        tandem::Repeat max_repeat {};
        Sequence motif(3, 'N');
        for (const auto& repeat : repeats) {
            set_motif(haplotype, repeat, motif);
            const auto penalty = get_open_penalty(motif, repeat.length);
            fill_n_if_less(std::next(std::begin(gap_open_penalities), repeat.pos), repeat.length, penalty);
            if (repeat.length > max_repeat.length) {
                max_repeat = repeat;
            }
        }
        set_motif(haplotype, max_repeat, motif);
        gap_extend_penalty = get_extension_penalty(motif, max_repeat.length);
    } else {
        gap_extend_penalty = get_default_extension_penalty();
    }
}
    
} // namespace octopus
