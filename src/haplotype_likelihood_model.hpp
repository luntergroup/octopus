//
//  haplotype_likelihood_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef haplotype_likelihood_model_hpp
#define haplotype_likelihood_model_hpp

#include <vector>
#include <iterator>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <functional>
#include <stdexcept>

#include <boost/optional.hpp>

#include "common.hpp"
#include "contig_region.hpp"
#include "haplotype.hpp"
#include "snv_error_model.hpp"
#include "indel_error_model.hpp"
#include "pair_hmm.hpp"

#include "timers.hpp"

class AlignedRead;

namespace Octopus
{
class HaplotypeLikelihoodModel
{
public:
    using PenaltyType = PairHMM::Model::PenaltyType;
    
    struct FlankState
    {
        ContigRegion::SizeType lhs_flank, rhs_flank;
    };
    
    class ShortHaplotypeError;
    
    using MappingPosition       = std::size_t;
    using MappingPositionVector = std::vector<MappingPosition>;
    using MappingPositionItr    = MappingPositionVector::const_iterator;
    
    HaplotypeLikelihoodModel();
    
    HaplotypeLikelihoodModel(SnvErrorModel snv_model, IndelErrorModel indel_model);
    
    HaplotypeLikelihoodModel(SnvErrorModel snv_model, IndelErrorModel indel_model,
                             const Haplotype& haplotype,
                             boost::optional<FlankState> flank_state = boost::none);
    
    HaplotypeLikelihoodModel(const HaplotypeLikelihoodModel&)            = default;
    HaplotypeLikelihoodModel& operator=(const HaplotypeLikelihoodModel&) = default;
    HaplotypeLikelihoodModel(HaplotypeLikelihoodModel&&)                 = default;
    HaplotypeLikelihoodModel& operator=(HaplotypeLikelihoodModel&&)      = default;
    
    ~HaplotypeLikelihoodModel() = default;
    
    void reset(const Haplotype& haplotype, boost::optional<FlankState> flank_state = boost::none);
    void clear() noexcept;
    
    // ln p(read | haplotype, model)
    
    double ln_probability(const AlignedRead& read) const;
    
    double ln_probability(const AlignedRead& read,
                          const MappingPositionVector& mapping_positions) const;
    
    double ln_probability(const AlignedRead& read,
                          MappingPositionItr first_mapping_position,
                          MappingPositionItr last_mapping_position) const;
    
private:
    SnvErrorModel snv_error_model_;
    IndelErrorModel indel_error_model_;
    
    const Haplotype* haplotype_;
    
    boost::optional<FlankState> haplotype_flank_state_;
    
    std::vector<char> haplotype_snv_forward_mask_, haplotype_snv_reverse_mask_;
    std::vector<PenaltyType> haplotype_snv_forward_priors_, haplotype_snv_reverse_priors_;
    
    std::vector<PenaltyType> haplotype_gap_open_penalities_;
    PenaltyType haplotype_gap_extension_penalty_;
};

class HaplotypeLikelihoodModel::ShortHaplotypeError : public std::runtime_error
{
public:
    using SizeType = Haplotype::SizeType;
    
    ShortHaplotypeError() = delete;
    
    ShortHaplotypeError(const Haplotype& haplotype, SizeType required_extension);
    
    const Haplotype& haplotype() const noexcept;
    
    SizeType required_extension() const noexcept;
    
private:
    const Haplotype& haplotype_;
    
    SizeType required_extension_;
};
} // namespace Octopus

#endif /* haplotype_likelihood_model_hpp */
