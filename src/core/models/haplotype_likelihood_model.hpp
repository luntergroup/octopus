// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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

#include <config/common.hpp>
#include <basics/contig_region.hpp>
#include <core/types/haplotype.hpp>

#include "readerror/snv_error_model.hpp"
#include "readerror/indel_error_model.hpp"
#include "pairhmm/pair_hmm.hpp"

#include "timers.hpp"

namespace octopus {

class AlignedRead;    

class HaplotypeLikelihoodModel
{
public:
    using Penalty = hmm::Model::Penalty;
    
    struct FlankState
    {
        ContigRegion::Position lhs_flank, rhs_flank;
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
    
    static unsigned pad_requirement() noexcept;
    
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
    std::vector<Penalty> haplotype_snv_forward_priors_, haplotype_snv_reverse_priors_;
    
    std::vector<Penalty> haplotype_gap_open_penalities_;
    Penalty haplotype_gap_extension_penalty_;
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

} // namespace octopus

#endif /* haplotype_likelihood_model_hpp */
