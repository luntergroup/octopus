// Copyright (c) 2017 Daniel Cooke
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
#include "core/types/haplotype.hpp"
#include "mutation/snv_error_model.hpp"
#include "mutation/indel_error_model.hpp"
#include "pairhmm/pair_hmm.hpp"

#include "timers.hpp"

namespace octopus {

class AlignedRead;    

class HaplotypeLikelihoodModel
{
public:
    using Penalty = hmm::MutationModel::Penalty;
    
    struct FlankState
    {
        ContigRegion::Position lhs_flank, rhs_flank;
    };
    
    class ShortHaplotypeError;
    
    using MappingPosition       = std::size_t;
    using MappingPositionVector = std::vector<MappingPosition>;
    using MappingPositionItr    = MappingPositionVector::const_iterator;
    
    HaplotypeLikelihoodModel();
    
    HaplotypeLikelihoodModel(bool use_mapping_quality);
    
    HaplotypeLikelihoodModel(std::unique_ptr<SnvErrorModel> snv_model,
                             std::unique_ptr<IndelErrorModel> indel_model,
                             bool use_mapping_quality = true);
    
    HaplotypeLikelihoodModel(std::unique_ptr<SnvErrorModel> snv_model,
                             std::unique_ptr<IndelErrorModel> indel_model,
                             const Haplotype& haplotype,
                             boost::optional<FlankState> flank_state = boost::none,
                             bool use_mapping_quality = true);
    
    HaplotypeLikelihoodModel(const HaplotypeLikelihoodModel&)            = delete;
    HaplotypeLikelihoodModel& operator=(const HaplotypeLikelihoodModel&) = delete;
    HaplotypeLikelihoodModel(HaplotypeLikelihoodModel&&)                 = default;
    HaplotypeLikelihoodModel& operator=(HaplotypeLikelihoodModel&&)      = default;
    
    ~HaplotypeLikelihoodModel() = default;
    
    static unsigned pad_requirement() noexcept;
    
    void reset(const Haplotype& haplotype, boost::optional<FlankState> flank_state = boost::none);
    
    void clear() noexcept;
    
    // ln p(read | haplotype, model)
    double evaluate(const AlignedRead& read) const;
    double evaluate(const AlignedRead& read, const MappingPositionVector& mapping_positions) const;
    double evaluate(const AlignedRead& read,
                    MappingPositionItr first_mapping_position,
                    MappingPositionItr last_mapping_position) const;
    
private:
    std::unique_ptr<SnvErrorModel> snv_error_model_;
    std::unique_ptr<IndelErrorModel> indel_error_model_;
    
    const Haplotype* haplotype_;
    
    boost::optional<FlankState> haplotype_flank_state_;
    
    std::vector<char> haplotype_snv_forward_mask_, haplotype_snv_reverse_mask_;
    std::vector<Penalty> haplotype_snv_forward_priors_, haplotype_snv_reverse_priors_;
    
    std::vector<Penalty> haplotype_gap_open_penalities_;
    Penalty haplotype_gap_extension_penalty_;
    bool use_mapping_quality_ = true;
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

HaplotypeLikelihoodModel make_haplotype_likelihood_model(const std::string sequencer, bool use_mapping_quality = true);

} // namespace octopus

#endif
