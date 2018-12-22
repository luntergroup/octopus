// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "custom_repeat_based_indel_error_model.hpp"

#include <utility>
#include <algorithm>
#include <iterator>
#include <array>

#include <boost/numeric/conversion/cast.hpp>
#include <boost/lexical_cast.hpp>

#include <iostream>

namespace octopus {

std::unique_ptr<IndelErrorModel> CustomRepeatBasedIndelErrorModel::do_clone() const
{
    return std::make_unique<CustomRepeatBasedIndelErrorModel>(*this);
}

namespace {

template <typename C, typename T>
static auto get_min_penalty(const C& penalties, const T length) noexcept
{
    return (length < penalties.size()) ? penalties[length] : penalties.back();
}

} // namespace

CustomRepeatBasedIndelErrorModel::CustomRepeatBasedIndelErrorModel(MotifPenaltyMap gap_open_penalties, PenaltyType extend_penalty)
: gap_open_penalties_ {std::move(gap_open_penalties)}
, gap_extend_penalties_ {}
, default_gap_open_ {get_min_penalty(gap_open_penalties_.at("A"), 0u)}
, default_gap_extend_ {extend_penalty}
, ns_ {}
{
    ns_.reserve(11);
    for (std::size_t i {0}; i <= 10; ++i) ns_.emplace_back(i, 'N');
}

CustomRepeatBasedIndelErrorModel::CustomRepeatBasedIndelErrorModel(MotifPenaltyMap gap_open_penalties, MotifPenaltyMap gap_extend_penalties)
: gap_open_penalties_ {std::move(gap_open_penalties)}
, gap_extend_penalties_ {std::move(gap_extend_penalties)}
, default_gap_open_ {get_min_penalty(gap_open_penalties_.at("A"), 0u)}
, default_gap_extend_ {get_min_penalty(gap_extend_penalties_->at("A"), 0u)}
, ns_ {}
{
    ns_.reserve(11);
    for (std::size_t i {0}; i <= 10; ++i) ns_.emplace_back(i, 'N');
}

CustomRepeatBasedIndelErrorModel::PenaltyType
CustomRepeatBasedIndelErrorModel::get_default_open_penalty() const noexcept
{
    return default_gap_open_;
}

CustomRepeatBasedIndelErrorModel::PenaltyType
CustomRepeatBasedIndelErrorModel::get_open_penalty(const Sequence& motif, const unsigned length) const noexcept
{
    const auto period = motif.size();
    auto itr = gap_open_penalties_.find(motif);
    if (itr == std::cend(gap_open_penalties_)) {
        itr = gap_open_penalties_.find(ns_[std::min(period, ns_.size() - 1)]);
    }
    return get_min_penalty(itr->second, length / period);
}

CustomRepeatBasedIndelErrorModel::PenaltyType
CustomRepeatBasedIndelErrorModel::get_default_extension_penalty() const noexcept
{
    return default_gap_extend_;
}

CustomRepeatBasedIndelErrorModel::PenaltyType
CustomRepeatBasedIndelErrorModel::get_extension_penalty(const Sequence& motif, const unsigned length) const noexcept
{
    if (!gap_extend_penalties_) return get_default_extension_penalty();
    const auto period = motif.size();
    auto itr = gap_extend_penalties_->find(motif);
    if (itr == std::cend(*gap_extend_penalties_)) {
        itr = gap_extend_penalties_->find(ns_[std::min(period, ns_.size() - 1)]);
    }
    return get_min_penalty(itr->second, length / period);
}

CustomRepeatBasedIndelErrorModel::MotifPenaltyMap make_penalty_map(std::string model)
{
    std::string token;
    CustomRepeatBasedIndelErrorModel::MotifPenaltyMap result {};
    for (auto token_itr = std::cbegin(model); token_itr != std::cend(model);) {
        auto motif_end_itr = std::find(token_itr, std::cend(model), ':');
        std::string motif {token_itr, motif_end_itr};
        token_itr = std::next(motif_end_itr);
        CustomRepeatBasedIndelErrorModel::PenaltyVector penalties {};
        penalties.reserve(100);
        for (; *std::prev(token_itr) != '\n'; ++token_itr) {
            constexpr std::array<char, 2> delims {',', '\n'};
            const auto token_end_itr = std::find_first_of(token_itr, std::cend(model), std::cbegin(delims), std::cend(delims));
            token.assign(token_itr, token_end_itr);
            try {
                penalties.push_back(boost::numeric_cast<CustomRepeatBasedIndelErrorModel::PenaltyType>(boost::lexical_cast<int>(token)));
            } catch (const boost::bad_lexical_cast&) {
                throw; // TODO: throw meaningful exception
            }
            token_itr = token_end_itr;
            if (token_itr == std::cend(model)) break;
        }
        result.emplace(std::move(motif), std::move(penalties));
    }
    return result;
}

} // namespace octopus