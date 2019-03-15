// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef haplotype_likelihood_model_hpp
#define haplotype_likelihood_model_hpp

#include <vector>
#include <iterator>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <memory>
#include <stdexcept>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/contig_region.hpp"
#include "basics/cigar_string.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/error/snv_error_model.hpp"
#include "core/models/error/indel_error_model.hpp"
#include "pairhmm/pair_hmm.hpp"

#include "timers.hpp"

namespace octopus {

class HaplotypeLikelihoodModel
{
public:
    using LogProbability = double;
    using Penalty = hmm::MutationModel::Penalty;
    
    struct Config
    {
        bool use_mapping_quality = true;
        boost::optional<AlignedRead::MappingQuality> mapping_quality_cap_trigger = boost::none;
        AlignedRead::MappingQuality mapping_quality_cap = 120;
        bool use_flank_state = true;
    };
    
    struct FlankState
    {
        ContigRegion::Position lhs_flank, rhs_flank;
    };
    
    class ShortHaplotypeError;
    
    using MappingPosition       = std::size_t;
    using MappingPositionVector = std::vector<MappingPosition>;
    using MappingPositionItr    = MappingPositionVector::const_iterator;
    
    struct Alignment
    {
        MappingPosition mapping_position;
        CigarString cigar;
        LogProbability likelihood;
    };
    
    HaplotypeLikelihoodModel();
    HaplotypeLikelihoodModel(Config config);
    HaplotypeLikelihoodModel(std::unique_ptr<SnvErrorModel> snv_model,
                             std::unique_ptr<IndelErrorModel> indel_model);
    HaplotypeLikelihoodModel(std::unique_ptr<SnvErrorModel> snv_model,
                             std::unique_ptr<IndelErrorModel> indel_model,
                             Config config);
    
    HaplotypeLikelihoodModel(const HaplotypeLikelihoodModel&);
    HaplotypeLikelihoodModel& operator=(const HaplotypeLikelihoodModel&);
    HaplotypeLikelihoodModel(HaplotypeLikelihoodModel&&)            = default;
    HaplotypeLikelihoodModel& operator=(HaplotypeLikelihoodModel&&) = default;
    
    friend void swap(HaplotypeLikelihoodModel& lhs, HaplotypeLikelihoodModel& rhs) noexcept;
    
    ~HaplotypeLikelihoodModel() = default;
    
    static unsigned pad_requirement() noexcept;
    
    bool can_use_flank_state() const noexcept;
    
    void reset(const Haplotype& haplotype, boost::optional<FlankState> flank_state = boost::none);
    
    void clear() noexcept;
    
    // ln p(read | haplotype, model)
    LogProbability evaluate(const AlignedRead& read) const;
    LogProbability evaluate(const AlignedRead& read, const MappingPositionVector& mapping_positions) const;
    LogProbability evaluate(const AlignedRead& read, MappingPositionItr first_mapping_position, MappingPositionItr last_mapping_position) const;
    
    Alignment align(const AlignedRead& read) const;
    Alignment align(const AlignedRead& read, const MappingPositionVector& mapping_positions) const;
    Alignment align(const AlignedRead& read, MappingPositionItr first_mapping_position, MappingPositionItr last_mapping_position) const;
    
private:
    std::unique_ptr<SnvErrorModel> snv_error_model_;
    std::unique_ptr<IndelErrorModel> indel_error_model_;
    
    const Haplotype* haplotype_;
    
    boost::optional<FlankState> haplotype_flank_state_;
    
    std::vector<char> haplotype_snv_forward_mask_, haplotype_snv_reverse_mask_;
    std::vector<Penalty> haplotype_snv_forward_priors_, haplotype_snv_reverse_priors_;
    
    std::vector<Penalty> haplotype_gap_open_penalities_, haplotype_gap_extend_penalities_;
    Config config_;
};

class HaplotypeLikelihoodModel::ShortHaplotypeError : public std::runtime_error
{
public:
    using Length = Haplotype::NucleotideSequence::size_type;
    
    ShortHaplotypeError() = delete;
    
    ShortHaplotypeError(const Haplotype& haplotype, Length required_extension);
    
    const Haplotype& haplotype() const noexcept;
    
    Length required_extension() const noexcept;
    
private:
    const Haplotype& haplotype_;
    Length required_extension_;
};

HaplotypeLikelihoodModel make_haplotype_likelihood_model(const std::string label, bool use_mapping_quality = true);

} // namespace octopus

#endif
