// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "hiseq_indel_error_model.hpp"

namespace octopus {

std::unique_ptr<IndelErrorModel> HiSeqIndelErrorModel::do_clone() const
{
    return std::make_unique<HiSeqIndelErrorModel>(*this);
}

namespace {

static constexpr std::array<IndelErrorModel::PenaltyType, 50> AT_homopolymer_penalities =
{{
 // 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
 60,60,45,43,41,38,35,32,29,25,21,19,18,18,17,16,16,15,15,14,
 14,13,12,11,11,10,9,9,8,7,7,7,6,6,6,6,6,6,6,5,4,4,4,4,3,3,3,3,3,3
 }};
static constexpr std::array<IndelErrorModel::PenaltyType, 50> CG_homopolymer_penalties =
{{
 // 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
 60,60,45,42,39,34,30,24,21,19,16,13,12,10,8,8,8,7,7,6,6,6,6,
 6,6,5,5,5,5,5,5,4,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3
 }};
static constexpr std::array<IndelErrorModel::PenaltyType, 50> dinucleotide_repeat_penalties =
{{
 // 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
 60,60,42,40,35,30,26,24,22,21,20,19,18,18,17,17,16,16,15,15,
 15,14,13,13,12,12,11,10,10,10,9,9,9,7,7,6,4,4,4,4,3,3,3,3,3,3,3,3,3,3
 }};
static constexpr std::array<IndelErrorModel::PenaltyType, 50> trinucleotide_repeat_penalties =
{{
 // 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
 60,60,40,36,30,28,26,25,23,23,23,22,21,21,20,19,18,17,17,16,
 15,15,15,14,13,12,11,11,10,8,7,6,5,5,5,5,5,4,4,3,3,3,3,3,3,3,3,3,3,3
 }};
static constexpr std::array<IndelErrorModel::PenaltyType, 50> polynucleotide_repeat_penalties =
{{
 // 1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19
 60,60,50,46,42,38,32,28,26,25,24,23,22,21,18,17,17,16,15,14,
 13,12,11,10,9,8,7,6,6,6,5,5,5,4,4,4,3,3,3,3,3,3,3,3,3,3,3,3,3,3
 }};

template <typename C, typename T>
static auto get_min_penalty(const C& penalties, const T length) noexcept
{
    return (length < penalties.size()) ? penalties[length] : penalties.back();
}

} // namespace

HiSeqIndelErrorModel::PenaltyType
HiSeqIndelErrorModel::get_default_open_penalty() const noexcept
{
    return AT_homopolymer_penalities.front();
}

HiSeqIndelErrorModel::PenaltyType
HiSeqIndelErrorModel::get_open_penalty(const Sequence& motif, const unsigned length) const noexcept
{
    const auto period = motif.size();
    const auto periodicity = length / period;
    switch (period) {
        case 1:
        {
            if (motif[0] == 'A' || motif[0] == 'T') {
                return get_min_penalty(AT_homopolymer_penalities, periodicity);
            } else {
                return get_min_penalty(CG_homopolymer_penalties, periodicity);
            }
        }
        case 2:
        {
            auto result = get_min_penalty(dinucleotide_repeat_penalties, periodicity);
            if (result > 7 && (motif == "CG" || motif == "GC")) result -= 2;
            return result;
        }
        case 3:
        {
            return get_min_penalty(trinucleotide_repeat_penalties, periodicity);
        }
        default:
            return get_min_penalty(polynucleotide_repeat_penalties, periodicity);
    }
}

HiSeqIndelErrorModel::PenaltyType
HiSeqIndelErrorModel::get_default_extension_penalty() const noexcept
{
    return 3;
}

HiSeqIndelErrorModel::PenaltyType
HiSeqIndelErrorModel::get_extension_penalty(const Sequence& motif, const unsigned length) const noexcept
{
    switch (motif.size()) {
        case 2:
        case 3: return 2;
        default: return get_default_extension_penalty();
    }
}

} // namespace octopus
