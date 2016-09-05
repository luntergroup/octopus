// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_model_hpp
#define denovo_model_hpp

#include "core/types/haplotype.hpp"

namespace octopus {

class DeNovoModel
{
public:
    struct Parameters
    {
        double mutation_rate;
    };
    
    DeNovoModel() = delete;
    
    DeNovoModel(Parameters params);
    
    DeNovoModel(const DeNovoModel&) = default;
    DeNovoModel& operator=(const DeNovoModel&) = default;
    DeNovoModel(DeNovoModel&&) = default;
    DeNovoModel& operator=(DeNovoModel&&)      = default;
    
    ~DeNovoModel() = default;
    
    double evaluate(const Haplotype& target, const Haplotype& given) const;

private:
    Parameters params_;
};

} // namespace octopus

 
#endif
