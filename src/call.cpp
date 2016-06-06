//
//  call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "call.hpp"

#include <algorithm>

namespace Octopus
{
double Call::quality() const noexcept
{
    return quality_;
}

Call::GenotypeCall& Call::get_genotype_call(const SampleIdType& sample)
{
    return genotype_calls_.at(sample);
}

const Call::GenotypeCall& Call::get_genotype_call(const SampleIdType& sample) const
{
    return genotype_calls_.at(sample);
}

bool Call::is_phased(const SampleIdType& sample) const
{
    return static_cast<bool>(genotype_calls_.at(sample).phase);
}

bool Call::all_phased() const noexcept
{
    return std::all_of(std::cbegin(genotype_calls_), std::cend(genotype_calls_),
                       [] (const auto& p) {
                           return static_cast<bool>(p.second.phase);
                       });
}

void Call::set_phase(const SampleIdType& sample, PhaseCall phase)
{
    genotype_calls_.at(sample).phase = std::move(phase);
}

void Call::replace(const char old_base, const char replacement_base)
{
    this->replace_called_alleles(old_base, replacement_base);
    
    for (auto& p : genotype_calls_) {
        auto& called_genotype = p.second.genotype;
        
        auto it = std::find_if_not(std::cbegin(called_genotype), std::cend(called_genotype),
                                   [old_base] (const Allele& allele) {
                                       const auto& seq = allele.sequence();
                                       return std::find(std::cbegin(seq), std::cend(seq),
                                                        old_base) == std::cend(seq);
                                   });
        
        if (it != std::cend(called_genotype)) {
            Genotype<Allele> new_genotype {called_genotype.ploidy()};
            
            std::for_each(std::cbegin(called_genotype), it,
                          [&new_genotype] (const Allele& allele) {
                              new_genotype.emplace(allele);
                          });
            
            std::for_each(it, std::cend(called_genotype),
                          [&new_genotype, old_base, replacement_base] (const Allele& allele) {
                              Allele::SequenceType new_sequence {};
                              
                              const auto& old_sequence = allele.sequence();
                              
                              new_sequence.reserve(old_sequence.size());
                              
                              std::replace_copy(std::cbegin(old_sequence), std::cend(old_sequence),
                                                std::back_inserter(new_sequence),
                                                old_base, replacement_base);
                              
                              new_genotype.emplace(Allele {allele.mapped_region(), std::move(new_sequence)});
                          });
            
            called_genotype = std::move(new_genotype);
        }
    }
}

void Call::set_dummy_model_bayes_factor(const double bf) noexcept
{
    dummy_model_bayes_factor_ = bf;
}

boost::optional<double> Call::dummy_model_bayes_factor() const noexcept
{
    return dummy_model_bayes_factor_;
}
} // namespace Octopus
