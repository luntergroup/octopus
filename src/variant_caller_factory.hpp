//
//  variant_caller_factory.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef variant_caller_factory_hpp
#define variant_caller_factory_hpp

#include <unordered_map>
#include <memory>

#include "genomic_region.hpp"
#include "variant_caller.hpp"
#include "variant_caller_builder.hpp"

class ReferenceGenome;
class ReadPipe;
class CandidateGeneratorBuilder;

namespace octopus
{
    class VariantCallerFactory
    {
    public:
        using ContigName = GenomicRegion::ContigName;
        
        VariantCallerFactory() = delete;
        
        VariantCallerFactory(VariantCallerBuilder template_builder, unsigned default_ploidy);
        
        ~VariantCallerFactory() = default;
        
        VariantCallerFactory(const VariantCallerFactory&)            = default;
        VariantCallerFactory& operator=(const VariantCallerFactory&) = default;
        VariantCallerFactory(VariantCallerFactory&&)                 = default;
        VariantCallerFactory& operator=(VariantCallerFactory&&)      = default;
        
        VariantCallerFactory& set_reference(const ReferenceGenome& reference) noexcept;
        VariantCallerFactory& set_read_pipe(ReadPipe& read_pipe) noexcept;
        VariantCallerFactory& set_candidate_generator_builder(const CandidateGeneratorBuilder& generator) noexcept;
        
        VariantCallerFactory& set_contig_ploidy(const ContigName& contig, unsigned ploidy);
        
        std::unique_ptr<VariantCaller> make(const ContigName& contig) const;
        
    private:
        mutable VariantCallerBuilder template_builder_;
        std::unordered_map<ContigName, unsigned> contig_ploidies_;
        unsigned default_ploidy_;
    };
} // namespace octopus

#endif /* variant_caller_factory_hpp */
