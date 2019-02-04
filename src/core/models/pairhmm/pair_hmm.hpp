// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef pair_hmm_hpp
#define pair_hmm_hpp

#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>

#include "basics/cigar_string.hpp"
#include "exceptions/program_error.hpp"

namespace octopus { namespace hmm {

// This is the minimum number of bases the truth must exceed the target either side around
// the mapped position
unsigned min_flank_pad() noexcept;

struct MutationModel
{
    using Penalty = std::int8_t;
    const std::vector<char>& snv_mask;
    const std::vector<Penalty>& snv_priors;
    const std::vector<Penalty>& gap_open;
    const std::vector<Penalty>& gap_extend;
    short nuc_prior = 2;
    std::size_t lhs_flank_size = 0, rhs_flank_size = 0;
};

struct VariableGapExtendMutationModel
{
    using Penalty = std::int8_t;
    using PenaltyVector = std::vector<Penalty>;
    Penalty mutation;
    const PenaltyVector& gap_open;
    const std::vector<Penalty>& gap_extend;
    short nuc_prior = 2;
};

struct VariableGapOpenMutationModel
{
    using Penalty = std::int8_t;
    using PenaltyVector = std::vector<Penalty>;
    Penalty mutation;
    const PenaltyVector& gap_open;
    short gap_extend;
    short nuc_prior = 2;
};

struct FlatGapMutationModel
{
    std::int8_t mutation;
    short gap_open, gap_extend;
    short nuc_prior = 2;
};

struct Alignment
{
    std::size_t target_offset;
    CigarString cigar;
    double likelihood;
};

class HMMOverflow : public ProgramError
{
public:
    using Sequence = std::string;
    
    HMMOverflow() = delete;
    HMMOverflow(const Sequence& target, const Sequence& truth) : target_ {target}, truth_ {truth} {}
    virtual ~HMMOverflow() override = default;
    
    const Sequence& target() const noexcept;
    const Sequence& truth() const noexcept;
    
private:
    const Sequence& target_, truth_;
    
    std::string do_why() const override { return "Pair HMM alignment overflowed"; }
    std::string do_where() const override { return "hmm::align"; }
};

// p(target | truth, target_qualities, target_offset, model)
//
// Warning: The target must be contained by the truth by at least
// min_flank_pad() on either side.
double evaluate(const std::string& target, const std::string& truth,
                const std::vector<std::uint8_t>& target_qualities,
                std::size_t target_offset,
                const MutationModel& model);

Alignment&
align(const std::string& target, const std::string& truth,
      const std::vector<std::uint8_t>& target_qualities,
      std::size_t target_offset,
      const MutationModel& model,
      Alignment& result);

Alignment
align(const std::string& target, const std::string& truth,
      const std::vector<std::uint8_t>& target_qualities,
      std::size_t target_offset,
      const MutationModel& model);

// p(target | truth, model)
//
// Warning: The target must be contained by the truth by exactly
// min_flank_pad() on either side.
double evaluate(const std::string& target, const std::string& truth,
                const VariableGapExtendMutationModel& model) noexcept;

Alignment&
align(const std::string& target, const std::string& truth,
      const VariableGapExtendMutationModel& model,
      Alignment& result);

Alignment
align(const std::string& target, const std::string& truth,
      const VariableGapExtendMutationModel& model);

// p(target | truth, model)
//
// Warning: The target must be contained by the truth by exactly
// min_flank_pad() on either side.
double evaluate(const std::string& target, const std::string& truth,
                const VariableGapOpenMutationModel& model) noexcept;

Alignment&
align(const std::string& target, const std::string& truth,
      const VariableGapOpenMutationModel& model,
      Alignment& result);

Alignment
align(const std::string& target, const std::string& truth,
      const VariableGapOpenMutationModel& model);

// p(target | truth, model)
//
// Warning: The target must be contained by the truth by exactly
// min_flank_pad() on either side.
double evaluate(const std::string& target, const std::string& truth,
                const FlatGapMutationModel& model) noexcept;

} // namespace hmm
} // namespace octopus

#endif
