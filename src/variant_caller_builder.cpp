//
//  variant_caller_builder.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_caller_builder.hpp"

#include "individual_caller.hpp"
#include "population_caller.hpp"
#include "cancer_caller.hpp"
#include "pedigree_caller.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
    // public methods
    
    VariantCallerBuilder::Parameters::Parameters(const ReferenceGenome& reference,
                                                 ReadPipe& read_pipe,
                                                 const CandidateGeneratorBuilder& candidate_generator_builder)
    :
    reference {reference},
    read_pipe {read_pipe},
    candidate_generator_builder {candidate_generator_builder}
    {}
    
    VariantCallerBuilder::VariantCallerBuilder(const ReferenceGenome& reference,
                                               ReadPipe& read_pipe,
                                               const CandidateGeneratorBuilder& candidate_generator_builder)
    :
    parameters_ {reference, read_pipe, candidate_generator_builder},
    factory_ {generate_factory()}
    {}
    
    VariantCallerBuilder::VariantCallerBuilder(const VariantCallerBuilder& other)
    :
    parameters_ {other.parameters_},
    factory_    {generate_factory()}
    {}
    
    VariantCallerBuilder& VariantCallerBuilder::operator=(const VariantCallerBuilder& other)
    {
        parameters_ = other.parameters_;
        factory_    = generate_factory();
        return *this;
    }
    
    VariantCallerBuilder::VariantCallerBuilder(VariantCallerBuilder&& other)
    :
    parameters_ {std::move(other.parameters_)},
    factory_    {generate_factory()}
    {}
    
    VariantCallerBuilder& VariantCallerBuilder::operator=(VariantCallerBuilder&& other)
    {
        std::swap(parameters_, other.parameters_);
        factory_ = generate_factory();
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_reference(const ReferenceGenome& reference) noexcept
    {
        parameters_.reference = reference;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_read_pipe(ReadPipe& read_pipe) noexcept
    {
        parameters_.read_pipe = read_pipe;
        return *this;
    }
    
    VariantCallerBuilder&
    VariantCallerBuilder::set_candidate_generator_builder(const CandidateGeneratorBuilder& candidate_generator_builder) noexcept
    {
        parameters_.candidate_generator_builder = candidate_generator_builder;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_ploidy(const unsigned ploidy) noexcept
    {
        parameters_.ploidy = ploidy;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_caller(std::string caller)
    {
        parameters_.caller = std::move(caller);
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_refcall_type(const VariantCaller::RefCallType refcall_type) noexcept
    {
        parameters_.refcall_type = refcall_type;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_sites_only() noexcept
    {
        parameters_.call_sites_only = true;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_min_variant_posterior(const double min_posterior) noexcept
    {
        parameters_.min_variant_posterior = min_posterior;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_min_refcall_posterior(const double min_posterior) noexcept
    {
        parameters_.min_refcall_posterior = min_posterior;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_max_haplotypes(const unsigned max_haplotypes) noexcept
    {
        parameters_.max_haplotypes = max_haplotypes;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_lagging(const bool allow_lagging) noexcept
    {
        parameters_.allow_lagging = allow_lagging;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_min_phase_score(const double min_phase_score) noexcept
    {
        parameters_.min_phase_score = min_phase_score;
        return *this;
    }
    
    // cancer
    VariantCallerBuilder& VariantCallerBuilder::set_normal_sample(SampleIdType normal_sample)
    {
        parameters_.normal_sample = std::move(normal_sample);
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_somatic_mutation_rate(double somatic_mutation_rate)
    {
        parameters_.somatic_mutation_rate = somatic_mutation_rate;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_min_somatic_posterior(const double min_posterior) noexcept
    {
        parameters_.min_somatic_posterior = min_posterior;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_somatic_only_calls() noexcept
    {
        parameters_.call_somatics_only = true;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_somatic_and_variant_calls() noexcept
    {
        parameters_.call_somatics_only = false;
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_somatic_and_variant_and_refcalls_calls() noexcept
    {
        parameters_.call_somatics_only = false;
        return *this;
    }
    
    // trio
    
    VariantCallerBuilder& VariantCallerBuilder::set_maternal_sample(SampleIdType mother)
    {
        parameters_.maternal_sample = std::move(mother);
        return *this;
    }
    
    VariantCallerBuilder& VariantCallerBuilder::set_paternal_sample(SampleIdType father)
    {
        parameters_.paternal_sample = std::move(father);
        return *this;
    }
    
    // pedigree
    
    VariantCallerBuilder& VariantCallerBuilder::set_pedigree(Pedigree pedigree)
    {
        parameters_.pedigree = std::move(pedigree);
        return *this;
    }
    
    // build
    
    std::unique_ptr<VariantCaller> VariantCallerBuilder::build() const
    {
        if (factory_.count(parameters_.caller) == 0) return nullptr;
        return factory_.at(parameters_.caller)();
    }
    
    // private methods
    
    VariantCallerBuilder::CallerFactoryMap VariantCallerBuilder::generate_factory() const
    {
        VariantCaller::CallerParameters general_parameters {
            parameters_.max_haplotypes,
            parameters_.refcall_type,
            parameters_.call_sites_only,
            parameters_.allow_lagging,
            parameters_.min_phase_score
        };
        
        return CallerFactoryMap {
            {"individual", [this, general_parameters = std::move(general_parameters)] () {
                return std::make_unique<IndividualVariantCaller>(parameters_.reference,
                                                                 parameters_.read_pipe,
                                                                 parameters_.candidate_generator_builder.get().build(),
                                                                 std::move(general_parameters),
                                                                 IndividualVariantCaller::CallerParameters {
                                                                     parameters_.min_variant_posterior,
                                                                     parameters_.min_refcall_posterior,
                                                                     parameters_.ploidy
                                                                 });
            }},
            {"population", [this, general_parameters = std::move(general_parameters)] () {
                return std::make_unique<PopulationVariantCaller>(parameters_.reference,
                                                                 parameters_.read_pipe,
                                                                 parameters_.candidate_generator_builder.get().build(),
                                                                 std::move(general_parameters),
                                                                 PopulationVariantCaller::CallerParameters {
                                                                     parameters_.min_variant_posterior,
                                                                     parameters_.min_refcall_posterior,
                                                                     parameters_.ploidy
                                                                 });
            }},
            {"cancer", [this, general_parameters = std::move(general_parameters)] () {
                return std::make_unique<CancerVariantCaller>(parameters_.reference,
                                                             parameters_.read_pipe,
                                                             parameters_.candidate_generator_builder.get().build(),
                                                             std::move(general_parameters),
                                                             CancerVariantCaller::CallerParameters {
                                                                 parameters_.min_variant_posterior,
                                                                 parameters_.min_somatic_posterior,
                                                                 parameters_.min_refcall_posterior,
                                                                 parameters_.ploidy,
                                                                 parameters_.normal_sample,
                                                                 parameters_.somatic_mutation_rate,
                                                                 parameters_.call_somatics_only
                                                             });
            }}
        };
    }
} // namespace Octopus
