// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef denovo_model_hpp
#define denovo_model_hpp

#include <unordered_map>

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
    DeNovoModel(Parameters parameters);
    
    DeNovoModel(const DeNovoModel&)            = default;
    DeNovoModel& operator=(const DeNovoModel&) = default;
    DeNovoModel(DeNovoModel&&)                 = default;
    DeNovoModel& operator=(DeNovoModel&&)      = default;
    
    ~DeNovoModel() = default;
    
    // ln p(target | given)
    double evaluate(const Haplotype& target, const Haplotype& given) const;
    
private:
    Parameters parameters_;
    
    mutable std::unordered_map<Haplotype, std::unordered_map<Haplotype, double>> cache_;
};

} // namespace octopus

 
#endif
