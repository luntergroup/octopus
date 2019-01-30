// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "basic_repeat_based_indel_error_model.hpp"

#include <vector>
#include <iterator>
#include <algorithm>

namespace octopus {

namespace {

template <typename T, std::size_t N>
void copy(const std::vector<T>& src, std::array<T, N>& dst) noexcept
{
    auto itr = std::copy(std::cbegin(src), std::next(std::cbegin(src), std::min(src.size(), N)), std::begin(dst));
    std::fill(itr, std::end(dst), src.back());
}

} // namespace

BasicRepeatBasedIndelErrorModel::BasicRepeatBasedIndelErrorModel(Parameters params)
{
    copy(params.AT_homopolymer_open_penalities, AT_homopolymer_open_penalities_);
    copy(params.CG_homopolymer_open_penalties, CG_homopolymer_open_penalties_);
    copy(params.dinucleotide_repeat_open_penalties, dinucleotide_repeat_open_penalties_);
    copy(params.trinucleotide_repeat_open_penalties, trinucleotide_repeat_open_penalties_);
    
    static const PenaltyVector homopolymer_extend_penalties {10, 10, 4, 4, 6, 6, 7, 8, 11, 13, 12, 11, 10, 9, 8, 7, 6, 6, 6, 5, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 3};
    static const PenaltyVector dinucleotide_extend_penalties {10, 10, 8, 5, 4, 2, 1};
    static const PenaltyVector trinucleotide_extend_penalties {10, 10, 6, 4, 3, 1};
    copy(homopolymer_extend_penalties, homopolymer_extend_penalties_);
    copy(dinucleotide_extend_penalties, dinucleotide_repeat_extend_penalties_);
    copy(trinucleotide_extend_penalties, trinucleotide_repeat_extend_penalties_);
    
    complex_open_penalty_ = dinucleotide_repeat_open_penalties_.front();
    complex_extend_penalty_ = 10;
}

std::unique_ptr<IndelErrorModel> BasicRepeatBasedIndelErrorModel::do_clone() const
{
    return std::make_unique<BasicRepeatBasedIndelErrorModel>(*this);
}

namespace {

template <typename C, typename T>
static auto get_min_penalty(const C& penalties, const T length) noexcept
{
    return (length < penalties.size()) ? penalties[length] : penalties.back();
}

} // namespace

BasicRepeatBasedIndelErrorModel::PenaltyType
BasicRepeatBasedIndelErrorModel::get_default_open_penalty() const noexcept
{
    return complex_open_penalty_;
}

BasicRepeatBasedIndelErrorModel::PenaltyType
BasicRepeatBasedIndelErrorModel::get_open_penalty(const Sequence& motif, const unsigned length) const noexcept
{
    const auto period = motif.size();
    const auto periodicity = length / period;
    switch (period) {
        case 1:
        {
            if (motif[0] == 'A' || motif[0] == 'T') {
                return get_min_penalty(AT_homopolymer_open_penalities_, periodicity);
            } else {
                return get_min_penalty(CG_homopolymer_open_penalties_, periodicity);
            }
        }
        case 2:
        {
            auto result = get_min_penalty(dinucleotide_repeat_open_penalties_, periodicity);
            if (result > 7 && (motif == "CG" || motif == "GC")) result -= 2;
            return result;
        }
        case 3:
        {
            return get_min_penalty(trinucleotide_repeat_open_penalties_, periodicity);
        }
        default:
            return get_min_penalty(trinucleotide_repeat_open_penalties_, periodicity);
    }
}

BasicRepeatBasedIndelErrorModel::PenaltyType
BasicRepeatBasedIndelErrorModel::get_default_extension_penalty() const noexcept
{
    return complex_extend_penalty_;
}

BasicRepeatBasedIndelErrorModel::PenaltyType
BasicRepeatBasedIndelErrorModel::get_extension_penalty(const Sequence& motif, const unsigned length) const noexcept
{
    const auto period = motif.size();
    const auto periodicity = length / period;
    switch (period) {
        case 1: return get_min_penalty(homopolymer_extend_penalties_, periodicity);
        case 2: return get_min_penalty(dinucleotide_repeat_extend_penalties_, periodicity);
        case 3:
        default: return get_min_penalty(trinucleotide_repeat_extend_penalties_, periodicity);
    }
}
    
} // namespace octopus
