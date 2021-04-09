// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indel_mutation_model_hpp
#define indel_mutation_model_hpp

#include <vector>
#include <cstdint>

#include "basics/tandem_repeat.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/variant.hpp"

namespace octopus {

class IndelMutationModel
{
public:
    using Probability = float;
    using ProbabilityVector = std::vector<Probability>;
    
    struct Parameters
    {
        double indel_mutation_prior;
        unsigned max_period = 10, max_periodicity = 50, max_indel_length = 50;
        double max_open_probability = 0.9, max_extend_probability = 1.0;
    };
    
    struct ContextIndelModel
    {
        ProbabilityVector gap_open;
        std::vector<ProbabilityVector> gap_extend;
    };
    
    IndelMutationModel() = delete;
    
    IndelMutationModel(Parameters params);
    
    IndelMutationModel(const IndelMutationModel&)            = default;
    IndelMutationModel& operator=(const IndelMutationModel&) = default;
    IndelMutationModel(IndelMutationModel&&)                 = default;
    IndelMutationModel& operator=(IndelMutationModel&&)      = default;
    
    ~IndelMutationModel() = default;
    
    ContextIndelModel evaluate(const Haplotype& haplotype) const;
    ContextIndelModel evaluate(const Haplotype& haplotype, std::vector<TandemRepeat>& repeats) const;
    
private:
    struct ModelCell
    {
        Probability open;
        ProbabilityVector extend;
    };
    using RepeatModel = std::vector<std::vector<ModelCell>>;
    
    Parameters params_;
    RepeatModel indel_repeat_model_;
};

IndelMutationModel::ContextIndelModel make_indel_model(const Haplotype& context, IndelMutationModel::Parameters params);
IndelMutationModel::ContextIndelModel make_indel_model(const Haplotype& context, IndelMutationModel::Parameters params, std::vector<TandemRepeat>& repeats);

IndelMutationModel::Probability
calculate_indel_probability(const IndelMutationModel::ContextIndelModel& model, std::size_t pos, std::size_t length) noexcept;

} // namespace octopus

#endif
