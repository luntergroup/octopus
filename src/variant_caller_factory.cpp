//
//  variant_caller_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 21/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_caller_factory.hpp"

#include <utility>

namespace Octopus
{
    VariantCallerFactory::VariantCallerFactory(VariantCallerBuilder template_builder,
                                               const unsigned default_ploidy)
    :
    template_builder_ {std::move(template_builder)},
    contig_ploidies_ {},
    default_ploidy_ {default_ploidy}
    {}
    
    void swap(VariantCallerFactory& lhs, VariantCallerFactory& rhs) noexcept
    {
        using std::swap;
        swap(lhs.template_builder_, rhs.template_builder_);
        swap(lhs.contig_ploidies_,  rhs.contig_ploidies_);
        swap(lhs.default_ploidy_,   rhs.default_ploidy_);
    }
    
    void VariantCallerFactory::set_reference(const ReferenceGenome& reference) noexcept
    {
        template_builder_.set_reference(reference);
    }
    
    void VariantCallerFactory::set_read_pipe(ReadPipe& read_pipe) noexcept
    {
        template_builder_.set_read_pipe(read_pipe);
    }
    
    void VariantCallerFactory::set_contig_ploidy(const ContigNameType& contig, const unsigned ploidy)
    {
        contig_ploidies_[contig] = ploidy;
    }
    
    std::unique_ptr<VariantCaller> VariantCallerFactory::make(const ContigNameType& contig) const
    {
        if (contig_ploidies_.count(contig) == 1) {
            template_builder_.set_ploidy(contig_ploidies_.at(contig));
        } else {
            template_builder_.set_ploidy(default_ploidy_);
        }
        return template_builder_.build();
    }
} // namespace Octopus
