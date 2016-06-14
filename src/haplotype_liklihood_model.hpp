//
//  haplotypeliklihood_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef haplotype_liklihood_model_hpp
#define haplotype_liklihood_model_hpp

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
#include "read_indel_error_model.hpp"
#include "pair_hmm.hpp"

#include "timers.hpp"

class AlignedRead;

namespace Octopus
{
class HaplotypeLikelihoodModel
{
public:
    using PenaltyType = ReadIndelErrorModel::PenaltyType;
    
    struct FlankState
    {
        ContigRegion::SizeType lhs_flank, rhs_flank;
    };
    
    class ShortHaplotypeError;
    
    using MapPositionItr = std::vector<std::size_t>::const_iterator;
    
    HaplotypeLikelihoodModel();
    
    HaplotypeLikelihoodModel(PenaltyType base_change_penalty, ReadIndelErrorModel indel_model);
    
    HaplotypeLikelihoodModel(PenaltyType base_change_penalty, ReadIndelErrorModel indel_model,
                             const Haplotype& haplotype, boost::optional<FlankState> flank_state);
    
    ~HaplotypeLikelihoodModel() = default;
    
    HaplotypeLikelihoodModel(const HaplotypeLikelihoodModel&)            = default;
    HaplotypeLikelihoodModel& operator=(const HaplotypeLikelihoodModel&) = default;
    HaplotypeLikelihoodModel(HaplotypeLikelihoodModel&&)                 = default;
    HaplotypeLikelihoodModel& operator=(HaplotypeLikelihoodModel&&)      = default;
    
    void set(const Haplotype& haplotype, boost::optional<FlankState> flank_state);
    
    void clear() noexcept;
    
    // ln p(read | haplotype, model)
    double log_probability(const AlignedRead& read,
                           MapPositionItr first_mapping_position,
                           MapPositionItr last_mapping_position) const;
    
private:
    PenaltyType base_change_penalty_;
    ReadIndelErrorModel indel_error_model_;
    
    const Haplotype* haplotype_;
    boost::optional<FlankState> haplotype_flank_state_;
    
    std::vector<std::int8_t> haplotype_gap_open_penalities_;
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

#endif /* haplotype_liklihood_model_hpp */
