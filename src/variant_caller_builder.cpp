//
//  variant_caller_builder.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_caller_builder.hpp"

#include <stdexcept>

#include "population_caller.hpp"
#include "cancer_caller.hpp"
#include "pedigree_caller.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
    // public methods
    
    VariantCallerBuilder::VariantCallerBuilder(const ReferenceGenome& reference,
                                               const CandidateGeneratorBuilder& candidate_generator_builder)
    :
    reference_ {reference},
    candidate_generator_builder_ {candidate_generator_builder},
    model_map_ {
        {"population", [&] () {
            return std::make_unique<PopulationVariantCaller>(reference_,
                                                             candidate_generator_builder_.get().build(),
                                                             refcall_type_,
                                                             min_variant_posterior_,
                                                             min_refcall_posterior_,
                                                             ploidy_);
        }},
        {"cancer", [&] () {
            return std::make_unique<CancerVariantCaller>(reference_,
                                                         candidate_generator_builder_.get().build(),
                                                         refcall_type_,
                                                         min_variant_posterior_,
                                                         min_somatic_posterior_,
                                                         min_refcall_posterior_,
                                                         normal_sample_.get(),
                                                         call_somatics_only_);
        }},
        {"trio", [&] () {
            return std::make_unique<PedigreeVariantCaller>(reference_,
                                                           candidate_generator_builder_.get().build(),
                                                           ploidy_,
                                                           maternal_sample_.get(),
                                                           paternal_sample_.get(),
                                                           min_variant_posterior_);
        }}
    }
    {}
    
    void VariantCallerBuilder::set_reference(const ReferenceGenome& reference)
    {
        reference_ = reference;
    }
    
    void VariantCallerBuilder::set_ploidy(const unsigned ploidy)
    {
        ploidy_ = ploidy;
    }
    
    void VariantCallerBuilder::set_model(std::string model)
    {
        model_ = std::move(model);
    }
    
    void VariantCallerBuilder::set_refcall_type(const VariantCaller::RefCallType refcall_type)
    {
        refcall_type_ = refcall_type;
    }
    
    void VariantCallerBuilder::set_min_variant_posterior(const double min_posterior)
    {
        min_variant_posterior_ = min_posterior;
    }
    
    void VariantCallerBuilder::set_min_refcall_posterior(const double min_posterior)
    {
        min_refcall_posterior_ = min_posterior;
    }
    
    // cancer
    void VariantCallerBuilder::set_normal_sample(SampleIdType normal_sample)
    {
        normal_sample_ = std::move(normal_sample);
    }
    
    void VariantCallerBuilder::set_min_somatic_posterior(const double min_posterior)
    {
        min_somatic_posterior_ = min_posterior;
    }
    
    void VariantCallerBuilder::set_somatic_only_calls()
    {
        call_somatics_only_ = true;
    }
    
    void VariantCallerBuilder::set_somatic_and_variant_calls()
    {
        call_somatics_only_ = false;
    }
    
    void VariantCallerBuilder::set_somatic_and_variant_and_refcalls_calls()
    {
        call_somatics_only_ = false;
    }
    
    // trio
    
    void VariantCallerBuilder::set_maternal_sample(SampleIdType mother)
    {
        maternal_sample_ = std::move(mother);
    }
    
    void VariantCallerBuilder::set_paternal_sample(SampleIdType father)
    {
        paternal_sample_ = std::move(father);
    }
    
    // pedigree
    
    void VariantCallerBuilder::set_pedigree(Pedigree pedigree)
    {
        pedigree_ = std::move(pedigree);
    }
    
    // build
    
    std::unique_ptr<VariantCaller> VariantCallerBuilder::build() const
    {
        if (model_map_.count(model_) == 0) {
            throw std::runtime_error {"VariantCallerBuilder: trying to build unknown model " + model_};
        }
        
        return model_map_.at(model_)();
    }
} // namespace Octopus
