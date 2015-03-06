//
//  variant_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_factory__
#define __Octopus__variant_factory__

#include "variant.h"
#include "utils.h"

class VariantFactory
{
public:
    using SizeType = Variant::SizeType;
    
    VariantFactory() = default;
    ~VariantFactory() = default;
    
    VariantFactory(const VariantFactory&)            = delete;
    VariantFactory& operator=(const VariantFactory&) = delete;
    VariantFactory(VariantFactory&&)                 = delete;
    VariantFactory& operator=(VariantFactory&&)      = delete;
    
    template <typename GenomicRegion_, typename StringType1, typename StringType2>
    Variant make(GenomicRegion_&& the_reference_allele_region, StringType1&& the_reference_allele,
                 StringType2&& the_alternative_allele) const;
    
    template <typename StringType1, typename StringType2, typename StringType3>
    Variant make(StringType1&& the_reference_contig_name, SizeType the_reference_begin,
                 StringType2&& the_reference_allele, StringType3&& the_alternative_allele) const;
    
private:
    template <typename T1, typename T2>
    std::function<double()> get_prior_probability_model(const T1& the_reference_allele,
                                                        const T2& the_alternative_allele) const;
};

template <typename GenomicRegion_, typename StringType1, typename StringType2>
Variant VariantFactory::make(GenomicRegion_&& the_reference_allele_region,
                             StringType1&& the_reference_allele,
                             StringType2&& the_alternative_allele) const
{
    auto&& prior_model = get_prior_probability_model(the_reference_allele, the_alternative_allele);
    return Variant {
        std::forward<GenomicRegion_>(the_reference_allele_region),
        std::forward<StringType1>(the_reference_allele),
        std::forward<StringType2>(the_alternative_allele),
        std::move(prior_model)
    };
}

template <typename StringType1, typename StringType2, typename StringType3>
Variant VariantFactory::make(StringType1&& the_reference_contig_name, SizeType the_reference_begin,
                             StringType2&& the_reference_allele, StringType3&& the_alternative_allele) const
{
    auto&& prior_model = get_prior_probability_model(the_reference_allele, the_alternative_allele);
    return Variant {
        std::forward<StringType1>(the_reference_contig_name),
        the_reference_begin,
        std::forward<StringType2>(the_reference_allele),
        std::forward<StringType3>(the_alternative_allele),
        std::move(prior_model)
    };
}

template <typename T1, typename T2>
std::function<double()>
VariantFactory::get_prior_probability_model(const T1& the_reference_allele,
                                            const T2& the_alternative_allele) const
{
    std::function<double()> prior_model {};
    if (stringlen(the_reference_allele) == stringlen(the_alternative_allele)) {
        if (stringlen(the_alternative_allele) == 1) {
            prior_model = [] () { return 1e-5; };
        } else {
            prior_model = [] () { return 1e-6; };
        }
    } else {
        if (stringlen(the_reference_allele) < stringlen(the_alternative_allele)) {
            prior_model = [] () { return 1e-7; };
        } else {
            prior_model = [] () { return 1e-8; };
        }
    }
    return prior_model;
}

#endif /* defined(__Octopus__variant_factory__) */
