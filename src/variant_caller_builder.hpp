//
//  variant_caller_builder.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef variant_caller_builder_hpp
#define variant_caller_builder_hpp

#include <unordered_map>
#include <string>
#include <memory>
#include <functional>

#include <boost/optional.hpp>

#include "common.hpp"
#include "variant_caller.hpp"
#include "read_pipe.hpp"
#include "candidate_generator_builder.hpp"
#include "pedigree.hpp"

namespace Octopus {
    
    class VariantCallerBuilder
    {
    public:
        VariantCallerBuilder()  = delete;
        
        explicit VariantCallerBuilder(const ReferenceGenome& reference,
                                      ReadPipe& read_pipe,
                                      const CandidateGeneratorBuilder& candidate_generator_builder);
        
        ~VariantCallerBuilder() = default;
        
        VariantCallerBuilder(const VariantCallerBuilder&);
        VariantCallerBuilder& operator=(const VariantCallerBuilder&);
        VariantCallerBuilder(VariantCallerBuilder&&);
        VariantCallerBuilder& operator=(VariantCallerBuilder&&);
        
        // common
        VariantCallerBuilder& set_reference(const ReferenceGenome& reference) noexcept;
        VariantCallerBuilder& set_read_pipe(ReadPipe& read_pipe) noexcept;
        VariantCallerBuilder& set_candidate_generator_builder(const CandidateGeneratorBuilder& candidate_generator_builder) noexcept;
        
        VariantCallerBuilder& set_ploidy(unsigned ploidy) noexcept;
        VariantCallerBuilder& set_model(std::string model);
        VariantCallerBuilder& set_refcall_type(VariantCaller::RefCallType refcall_type) noexcept;
        VariantCallerBuilder& set_min_variant_posterior(double min_posterior) noexcept;
        VariantCallerBuilder& set_min_refcall_posterior(double min_posterior) noexcept;
        
        // cancer
        VariantCallerBuilder& set_normal_sample(SampleIdType normal_sample);
        VariantCallerBuilder& set_min_somatic_posterior(double min_posterior) noexcept;
        VariantCallerBuilder& set_somatic_only_calls() noexcept;
        VariantCallerBuilder& set_somatic_and_variant_calls() noexcept;
        VariantCallerBuilder& set_somatic_and_variant_and_refcalls_calls() noexcept;
        
        // trio
        
        VariantCallerBuilder& set_maternal_sample(SampleIdType mother);
        VariantCallerBuilder& set_paternal_sample(SampleIdType father);
        
        // pedigree
        
        VariantCallerBuilder& set_pedigree(Pedigree pedigree);
        
        // build
        
        std::unique_ptr<VariantCaller> build() const;
        
    private:
        // common parameters
        std::reference_wrapper<const ReferenceGenome> reference_;
        std::reference_wrapper<ReadPipe> read_pipe_;
        
        unsigned ploidy_ = 2;
        std::string model_ = "population";
        std::reference_wrapper<const CandidateGeneratorBuilder> candidate_generator_builder_;
        VariantCaller::RefCallType refcall_type_ = VariantCaller::RefCallType::None;
        double min_variant_posterior_ = 0.99;
        double min_refcall_posterior_ = 0.99;
        
        // cancer
        
        boost::optional<SampleIdType> normal_sample_;
        double min_somatic_posterior_ = 0.99;
        bool call_somatics_only_ = true;
        
        // trio
        
        boost::optional<SampleIdType> maternal_sample_, paternal_sample_;
        
        // pedigree
        
        boost::optional<Pedigree> pedigree_;
        
        // factory
        
        using ModelFactoryMap = std::unordered_map<std::string, std::function<std::unique_ptr<VariantCaller>()>>;
        
        ModelFactoryMap factory_;
        
        ModelFactoryMap generate_factory() const;
    };
} // namespace Octopus

#endif /* variant_caller_builder_hpp */
