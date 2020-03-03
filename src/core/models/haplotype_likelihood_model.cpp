// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "haplotype_likelihood_model.hpp"

#include <utility>
#include <cmath>
#include <limits>
#include <cassert>

#include "core/models/error/error_model_factory.hpp"
#include "concepts/mappable.hpp"
#include "utils/maths.hpp"

namespace octopus {

HaplotypeLikelihoodModel::ShortHaplotypeError::ShortHaplotypeError(const Haplotype& haplotype,
                                                                   Length required_extension)
: std::runtime_error {"Haplotype is too short for alignment"}
, haplotype_ {haplotype}
, required_extension_ {required_extension}
{}

const Haplotype& HaplotypeLikelihoodModel::ShortHaplotypeError::haplotype() const noexcept
{
    return haplotype_;
}

HaplotypeLikelihoodModel::ShortHaplotypeError::Length
HaplotypeLikelihoodModel::ShortHaplotypeError::required_extension() const noexcept
{
    return required_extension_;
}

// public methods

const HaplotypeLikelihoodModel::Config& HaplotypeLikelihoodModel::config() const noexcept
{
    return config_;
}

void HaplotypeLikelihoodModel::set(Config config)
{
    if (config.max_indel_error != this->config_.max_indel_error) {
        if (config.use_int_scores) {
            hmm_ = HMM {config.max_indel_error, HMM::ScoreType::int32};
        } else {
            hmm_ = HMM {config.max_indel_error};
        }
    }
    config_ = std::move(config);
    if (config_.mapping_quality_cap_trigger && *config_.mapping_quality_cap_trigger >= config_.mapping_quality_cap) {
        config_.mapping_quality_cap_trigger = boost::none;
    }
}

unsigned HaplotypeLikelihoodModel::pad_requirement() const noexcept
{
    return hmm_.band_size();
}

void HaplotypeLikelihoodModel::reset(const Haplotype& haplotype, boost::optional<FlankState> flank_state)
{
    haplotype_ = std::addressof(haplotype);
    haplotype_flank_state_ = std::move(flank_state);
    if (snv_error_model_) {
        snv_error_model_->evaluate(haplotype,
                                   haplotype_snv_forward_mask_, haplotype_snv_forward_priors_,
                                   haplotype_snv_reverse_mask_, haplotype_snv_reverse_priors_);
    } else {
        // TODO: refactor HaplotypeLikelihoodModel to use another HMM evaluate overload without SNV model
        haplotype_snv_forward_priors_.assign(sequence_size(haplotype), 100);
        haplotype_snv_forward_mask_.assign(std::cbegin(haplotype.sequence()), std::cend(haplotype.sequence()));
        haplotype_snv_reverse_priors_.assign(sequence_size(haplotype), 100);
        haplotype_snv_reverse_mask_.assign(std::cbegin(haplotype.sequence()), std::cend(haplotype.sequence()));
    }
    if (indel_error_model_) {
        indel_error_model_->set_penalties(haplotype, haplotype_gap_open_penalities_, haplotype_gap_extend_penalities_);
    }
}

void HaplotypeLikelihoodModel::clear() noexcept
{
    haplotype_ = nullptr;
    haplotype_flank_state_ = boost::none;
}

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel()
: HaplotypeLikelihoodModel {make_snv_error_model(), make_indel_error_model(), Config {}} {}

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(Config config)
: HaplotypeLikelihoodModel {make_snv_error_model(), make_indel_error_model(), config} {}

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(std::unique_ptr<SnvErrorModel> snv_model,
                                                   std::unique_ptr<IndelErrorModel> indel_model)
: HaplotypeLikelihoodModel {std::move(snv_model), std::move(indel_model), Config {}} {}

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(std::unique_ptr<SnvErrorModel> snv_model,
                                                   std::unique_ptr<IndelErrorModel> indel_model,
                                                   Config config)
: snv_error_model_ {std::move(snv_model)}
, indel_error_model_ {std::move(indel_model)}
, haplotype_ {nullptr}
, haplotype_flank_state_ {}
, haplotype_gap_open_penalities_ {}
, haplotype_gap_extend_penalities_ {}
, config_ {config}
{
    if (config.use_int_scores) {
        hmm_ = HMM {config.max_indel_error, HMM::ScoreType::int32};
    } else {
        hmm_ = HMM {config.max_indel_error};
    }
    if (config_.mapping_quality_cap_trigger && *config_.mapping_quality_cap_trigger >= config_.mapping_quality_cap) {
        config_.mapping_quality_cap_trigger = boost::none;
    }
}

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(const HaplotypeLikelihoodModel& other)
{
    if (other.indel_error_model_) {
        indel_error_model_ = other.indel_error_model_->clone();
    } else {
        indel_error_model_ = nullptr;
    }
    if (other.snv_error_model_) {
        snv_error_model_ = other.snv_error_model_->clone();
    } else {
        snv_error_model_ = nullptr;
    }
    haplotype_ = other.haplotype_;
    haplotype_flank_state_ = other.haplotype_flank_state_;
    haplotype_snv_forward_mask_ = other.haplotype_snv_forward_mask_;
    haplotype_snv_reverse_mask_ = other.haplotype_snv_reverse_mask_;
    haplotype_snv_forward_priors_ = other.haplotype_snv_forward_priors_;
    haplotype_snv_reverse_priors_ = other.haplotype_snv_reverse_priors_;
    haplotype_gap_open_penalities_ = other.haplotype_gap_open_penalities_;
    haplotype_gap_extend_penalities_ = other.haplotype_gap_extend_penalities_;
    config_ = other.config_;
    hmm_ = other.hmm_;
}

HaplotypeLikelihoodModel& HaplotypeLikelihoodModel::operator=(const HaplotypeLikelihoodModel& other)
{
    if (this != &other) {
        HaplotypeLikelihoodModel tmp {other};
        std::swap(*this, tmp);
    }
    return *this;
}

void swap(HaplotypeLikelihoodModel& lhs, HaplotypeLikelihoodModel& rhs) noexcept
{
    using std::swap;
    swap(lhs.indel_error_model_, rhs.indel_error_model_);
    swap(lhs.snv_error_model_, rhs.snv_error_model_);
    swap(lhs.haplotype_, rhs.haplotype_);
    swap(lhs.haplotype_flank_state_, rhs.haplotype_flank_state_);
    swap(lhs.haplotype_snv_forward_mask_, rhs.haplotype_snv_forward_mask_);
    swap(lhs.haplotype_snv_reverse_mask_, rhs.haplotype_snv_reverse_mask_);
    swap(lhs.haplotype_snv_forward_priors_, rhs.haplotype_snv_forward_priors_);
    swap(lhs.haplotype_snv_reverse_priors_, rhs.haplotype_snv_reverse_priors_);
    swap(lhs.haplotype_gap_open_penalities_, rhs.haplotype_gap_open_penalities_);
    swap(lhs.haplotype_gap_extend_penalities_, rhs.haplotype_gap_extend_penalities_);
    swap(lhs.config_, rhs.config_);
    swap(lhs.hmm_, rhs.hmm_);
}

bool HaplotypeLikelihoodModel::can_use_flank_state() const noexcept
{
    return config_.use_flank_state;
}

HaplotypeLikelihoodModel::LogProbability
HaplotypeLikelihoodModel::evaluate(const AlignedRead& read) const
{
    const static MappingPositionVector empty {};
    return this->evaluate(read, std::cbegin(empty), std::cend(empty));
}

HaplotypeLikelihoodModel::LogProbability
HaplotypeLikelihoodModel::evaluate(const AlignedRead& read, const MappingPositionVector& mapping_positions) const
{
    return this->evaluate(read, std::cbegin(mapping_positions), std::cend(mapping_positions));
}

namespace {

template <typename pHMM>
int num_out_of_range_bases(const std::size_t mapping_position, const AlignedRead& read, const Haplotype& haplotype,
                           const pHMM& hmm) noexcept
{
    const auto required_pad = min_flank_pad(hmm);
    if (mapping_position < required_pad) {
        return required_pad - mapping_position;
    }
    const auto mapping_end = mapping_position + sequence_size(read) + required_pad;
    if (mapping_end > sequence_size(haplotype)) {
        return static_cast<int>(sequence_size(haplotype)) - static_cast<int>(mapping_end);
    } else {
        return 0;
    }
}

template <typename pHMM>
bool is_in_range(const std::size_t mapping_position, const AlignedRead& read, const Haplotype& haplotype, const pHMM& hmm) noexcept
{
    return num_out_of_range_bases(mapping_position, read, haplotype, hmm) == 0;
}

} // namespace

template <typename InputIt, typename pHMM>
HaplotypeLikelihoodModel::LogProbability
max_score(const AlignedRead& read, const Haplotype& haplotype,
          InputIt first_mapping_position, InputIt last_mapping_position,
          const pHMM& hmm)
{
    assert(contains(haplotype, read));
    using LogProbability = HaplotypeLikelihoodModel::LogProbability;
    using PositionType = typename std::iterator_traits<InputIt>::value_type;
    const auto original_mapping_position = static_cast<PositionType>(begin_distance(haplotype, read));
    auto max_log_probability = std::numeric_limits<LogProbability>::lowest();
    bool is_original_position_mapped {false}, has_in_range_mapping_position {false};
    std::for_each(first_mapping_position, last_mapping_position, [&] (const auto position) {
        if (position == original_mapping_position) {
            is_original_position_mapped = true;
        }
        if (is_in_range(position, read, haplotype, hmm)) {
            has_in_range_mapping_position = true;
            auto p = hmm.evaluate(read.sequence(), haplotype.sequence(), read.base_qualities(), position);
            max_log_probability = std::max(static_cast<LogProbability>(p), max_log_probability);
        }
    });
    if (!is_original_position_mapped && is_in_range(original_mapping_position, read, haplotype, hmm)) {
        has_in_range_mapping_position = true;
        auto p = hmm.evaluate(read.sequence(), haplotype.sequence(), read.base_qualities(), original_mapping_position);
        max_log_probability = std::max(static_cast<LogProbability>(p), max_log_probability);
    }
    if (!has_in_range_mapping_position) {
        const auto min_shift = num_out_of_range_bases(original_mapping_position, read, haplotype, hmm);
        auto final_mapping_position = original_mapping_position;
        if (min_shift > 0) {
            final_mapping_position += min_shift;
            if (!is_in_range(final_mapping_position, read, haplotype, hmm)) {
                throw HaplotypeLikelihoodModel::ShortHaplotypeError {haplotype, static_cast<unsigned>(min_shift)};
            }
        } else {
            const auto min_left_shift = static_cast<unsigned>(-min_shift);
            if (original_mapping_position >= min_left_shift) {
                final_mapping_position -= min_left_shift;
            } else {
                auto required_extension = min_left_shift - original_mapping_position;
                throw HaplotypeLikelihoodModel::ShortHaplotypeError {haplotype, required_extension};
            }
        }
        max_log_probability = hmm.evaluate(read.sequence(), haplotype.sequence(), read.base_qualities(), final_mapping_position);
    }
    assert(max_log_probability > std::numeric_limits<LogProbability>::lowest() && max_log_probability <= 0);
    return max_log_probability;
}

HaplotypeLikelihoodModel::LogProbability
HaplotypeLikelihoodModel::evaluate(const AlignedRead& read,
                                   MappingPositionItr first_mapping_position,
                                   MappingPositionItr last_mapping_position) const
{
    if (haplotype_ == nullptr) {
        throw std::runtime_error {"HaplotypeLikelihoodModel: no buffered Haplotype"};
    }
    const auto is_forward = !read.is_marked_reverse_mapped();
    HMM::ParameterType model {
        haplotype_gap_open_penalities_,
        haplotype_gap_extend_penalities_,
        is_forward ? haplotype_snv_forward_mask_ : haplotype_snv_reverse_mask_,
        is_forward ? haplotype_snv_forward_priors_ : haplotype_snv_reverse_priors_
    };
    if (haplotype_flank_state_) {
        model.lhs_flank_size = haplotype_flank_state_->lhs_flank;
        model.rhs_flank_size = haplotype_flank_state_->rhs_flank;
    } else {
        model.lhs_flank_size = 0;
        model.rhs_flank_size = 0;
    }
    hmm_.set(model);
    const auto ln_prob_given_mapped = max_score(read, *haplotype_, first_mapping_position, last_mapping_position, hmm_);
    if (config_.use_mapping_quality) {
        // This calculation is approximately
        // p(read | hap) = p(read missmapped) p(read | hap, missmapped)
        //                  + p(read correctly mapped) p(read | hap, correctly mapped)
        // = p(read correctly mapped) p(read | hap, correctly mapped)
        //      + p(read missmapped)
        // assuming p(read | hap, missmapped) = 1
        auto mapping_quality = read.mapping_quality();
        if (config_.mapping_quality_cap_trigger && mapping_quality >= *config_.mapping_quality_cap_trigger) {
            mapping_quality = config_.mapping_quality_cap;
        }
        using octopus::maths::constants::ln10Div10;
        const auto ln_prob_missmapped = -ln10Div10<> * mapping_quality;
        const auto ln_prob_mapped = std::log(1.0 - std::exp(ln_prob_missmapped));
        const auto result = maths::log_sum_exp(ln_prob_mapped + ln_prob_given_mapped, ln_prob_missmapped);
        return result > -1e-15 ? 0.0 : result;
    } else {
        return ln_prob_given_mapped  > -1e-15 ? 0.0 : ln_prob_given_mapped;
    }
}

HaplotypeLikelihoodModel::LogProbability
HaplotypeLikelihoodModel::evaluate(const AlignedTemplate& reads) const
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), LogProbability {0},
                           [this] (LogProbability curr, const AlignedRead& read) { return curr + this->evaluate(read); });
}

HaplotypeLikelihoodModel::LogProbability
HaplotypeLikelihoodModel::evaluate(const AlignedTemplate& reads, const std::vector<MappingPositionVector>& mapping_positions) const
{
    assert(reads.size() == mapping_positions.size());
    return std::inner_product(std::cbegin(reads), std::cend(reads), std::cbegin(mapping_positions), LogProbability {0},
                              std::plus<> {}, [this] (const AlignedRead& read, const MappingPositionVector& mapping_positions)
                              { return this->evaluate(read, mapping_positions); });
}

HaplotypeLikelihoodModel::Alignment
HaplotypeLikelihoodModel::align(const AlignedRead& read) const
{
    const static MappingPositionVector empty {};
    return this->align(read, std::cbegin(empty), std::cend(empty));
}

HaplotypeLikelihoodModel::Alignment
HaplotypeLikelihoodModel::align(const AlignedRead& read, const MappingPositionVector& mapping_positions) const
{
    return this->align(read, std::cbegin(mapping_positions), std::cend(mapping_positions));
}

template <typename InputIt, typename pHMM>
HaplotypeLikelihoodModel::Alignment
compute_optimal_alignment(const AlignedRead& read, const Haplotype& haplotype,
                          InputIt first_mapping_position, InputIt last_mapping_position,
                          const pHMM& hmm)
{
    assert(contains(haplotype, read));
    using LogProbability = HaplotypeLikelihoodModel::LogProbability;
    using PositionType = typename std::iterator_traits<InputIt>::value_type;
    const auto original_mapping_position = static_cast<PositionType>(begin_distance(haplotype, read));
    HaplotypeLikelihoodModel::Alignment result {};
    result.likelihood = std::numeric_limits<LogProbability>::lowest();
    bool is_original_position_mapped {false}, has_in_range_mapping_position {false};
    std::for_each(first_mapping_position, last_mapping_position, [&] (const auto position) {
        if (position == original_mapping_position) {
            is_original_position_mapped = true;
        }
        if (is_in_range(position, read, haplotype, hmm)) {
            has_in_range_mapping_position = true;
            auto alignment = hmm.align(read.sequence(), haplotype.sequence(), read.base_qualities(), position);
            if (alignment.likelihood > result.likelihood) {
                result.mapping_position = alignment.target_offset;
                result.likelihood = alignment.likelihood;
                result.cigar = std::move(alignment.cigar);
            }
        }
    });
    if (!is_original_position_mapped && is_in_range(original_mapping_position, read, haplotype, hmm)) {
        has_in_range_mapping_position = true;
        auto alignment = hmm.align(read.sequence(), haplotype.sequence(), read.base_qualities(), original_mapping_position);
        if (alignment.likelihood >= result.likelihood) {
            result.mapping_position = alignment.target_offset;
            result.likelihood = alignment.likelihood;
            result.cigar = std::move(alignment.cigar);
        }
    }
    if (!has_in_range_mapping_position) {
        const auto min_shift = num_out_of_range_bases(original_mapping_position, read, haplotype, hmm);
        auto final_mapping_position = original_mapping_position;
        if (min_shift > 0) {
            final_mapping_position += min_shift;
            if (!is_in_range(final_mapping_position, read, haplotype, hmm)) {
                throw HaplotypeLikelihoodModel::ShortHaplotypeError {haplotype, static_cast<unsigned>(min_shift)};
            }
        } else {
            const auto min_left_shift = static_cast<unsigned>(-min_shift);
            if (original_mapping_position >= min_left_shift) {
                final_mapping_position -= min_left_shift;
            } else {
                auto required_extension = min_left_shift - original_mapping_position;
                throw HaplotypeLikelihoodModel::ShortHaplotypeError {haplotype, required_extension};
            }
        }
        auto alignment = hmm.align(read.sequence(), haplotype.sequence(), read.base_qualities(), final_mapping_position);
        result.likelihood = alignment.likelihood;
        result.cigar = std::move(alignment.cigar);
        result.mapping_position = alignment.target_offset;
    }
    assert(result.likelihood > std::numeric_limits<LogProbability>::lowest() && result.likelihood <= 0);
    return result;
}

HaplotypeLikelihoodModel::Alignment
HaplotypeLikelihoodModel::align(const AlignedRead& read, MappingPositionItr first_mapping_position, MappingPositionItr last_mapping_position) const
{
    if (haplotype_ == nullptr) {
        throw std::runtime_error {"HaplotypeLikelihoodModel: no buffered Haplotype"};
    }
    const auto is_forward = !read.is_marked_reverse_mapped();
    HMM::ParameterType model {haplotype_gap_open_penalities_,
                              haplotype_gap_extend_penalities_,
                              is_forward ? haplotype_snv_forward_mask_ : haplotype_snv_reverse_mask_,
                              is_forward ? haplotype_snv_forward_priors_ : haplotype_snv_reverse_priors_};
    if (haplotype_flank_state_) {
        model.lhs_flank_size = haplotype_flank_state_->lhs_flank;
        model.rhs_flank_size = haplotype_flank_state_->rhs_flank;
    } else {
        model.lhs_flank_size = 0;
        model.rhs_flank_size = 0;
    }
    hmm_.set(model);
    auto result = compute_optimal_alignment(read, *haplotype_, first_mapping_position, last_mapping_position, hmm_);
    if (config_.use_mapping_quality) {
        auto mapping_quality = read.mapping_quality();
        if (config_.mapping_quality_cap_trigger && mapping_quality >= *config_.mapping_quality_cap_trigger) {
            mapping_quality = config_.mapping_quality_cap;
        }
        using octopus::maths::constants::ln10Div10;
        const auto ln_prob_missmapped = -ln10Div10<> * mapping_quality;
        const auto ln_prob_mapped = std::log(1.0 - std::exp(ln_prob_missmapped));
        result.likelihood = maths::log_sum_exp(ln_prob_mapped + result.likelihood, ln_prob_missmapped);
        result.likelihood = result.likelihood > -1e-15 ? 0.0 : result.likelihood;
    } else {
        result.likelihood = result.likelihood > -1e-15 ? 0.0 : result.likelihood;
    }
    return result;
}

HaplotypeLikelihoodModel make_haplotype_likelihood_model(const std::string label, bool use_mapping_quality)
{
    HaplotypeLikelihoodModel::Config config {};
    config.use_mapping_quality = use_mapping_quality;
    auto error_model = make_error_model(label);
    return HaplotypeLikelihoodModel {std::move(error_model.snv), std::move(error_model.indel), config};
}

} // namespace octopus
