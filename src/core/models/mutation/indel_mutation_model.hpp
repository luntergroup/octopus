// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indel_mutation_model_hpp
#define indel_mutation_model_hpp

#include <vector>
#include <cstdint>

#include "core/types/haplotype.hpp"
#include "core/types/variant.hpp"

namespace octopus {

class IndelMutationModel
{
public:
    struct Parameters
    {
        double indel_mutation_rate;
        unsigned max_period = 10, max_periodicity = 20;
        double max_open_rate = 1.0, max_extend_rate = 1.0;
    };
    
    struct ContextIndelModel
    {
        using Rate = double;
        using RateVector = std::vector<Rate>;
        RateVector gap_open, gap_extend;
    };
    
    IndelMutationModel() = delete;
    
    IndelMutationModel(Parameters params);
    
    IndelMutationModel(const IndelMutationModel&)            = default;
    IndelMutationModel& operator=(const IndelMutationModel&) = default;
    IndelMutationModel(IndelMutationModel&&)                 = default;
    IndelMutationModel& operator=(IndelMutationModel&&)      = default;
    
    ~IndelMutationModel() = default;
    
    ContextIndelModel evaluate(const Haplotype& haplotype) const;
    
private:
    struct ModelCell { double open, extend; };
    using RepeatModel = std::vector<std::vector<ModelCell>>;
    
    Parameters params_;
    RepeatModel indel_repeat_model_;
};

IndelMutationModel::ContextIndelModel make_indel_model(const Haplotype& context, IndelMutationModel::Parameters params);

} // namespace octopus

#endif
