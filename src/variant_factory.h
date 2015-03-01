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

class VariantFactory
{
public:
    using SizeType = Variant::SizeType;
    
    VariantFactory() = default;;
    
    VariantFactory(const VariantFactory&)            = delete;
    VariantFactory& operator=(const VariantFactory&) = delete;
    VariantFactory(VariantFactory&&)                 = delete;
    VariantFactory& operator=(VariantFactory&&)      = delete;
    
    template <typename T1, typename T2, typename T3>
    Variant make(T1&& the_region, T2&& sequence_removed, T3&& sequence_added) const;
    
    template <typename T1, typename T2, typename T3>
    Variant make(T1&& contig_name, SizeType contig_begin_pos, T2&& sequence_removed,
                 T3&& sequence_added) const;
};

template <typename T1, typename T2, typename T3>
Variant VariantFactory::make(T1&& the_region, T2&& sequence_removed, T3&& sequence_added) const
{
    std::function<double()> prior_model {};
    if (sequence_added.size() == sequence_removed.size()) {
        if (sequence_added.length() == 1) {
            prior_model = [] () { return 1e-5; };
        } else {
            prior_model = [] () { return 1e-6; };
        }
    } else {
        if (sequence_added.length() < sequence_removed.length()) {
            prior_model = [] () { return 1e-7; };
        } else {
            prior_model = [] () { return 1e-8; };
        }
    }
    return Variant {
        std::forward<T1>(the_region),
        std::forward<T2>(sequence_removed),
        std::forward<T3>(sequence_added),
        prior_model
    };
}

template <typename T1, typename T2, typename T3>
Variant VariantFactory::make(T1&& contig_name, SizeType contig_begin_pos,
                             T2&& sequence_removed, T3&& sequence_added) const
{
    std::function<double()> prior_model {};
    if (sequence_added.size() == sequence_removed.size()) {
        if (sequence_added.length() == 1) {
            prior_model = [] () { return 1e-5; };
        } else {
            prior_model = [] () { return 1e-6; };
        }
    } else {
        if (sequence_added.length() < sequence_removed.length()) {
            prior_model = [] () { return 1e-7; };
        } else {
            prior_model = [] () { return 1e-8; };
        }
    }
    return Variant {
        std::forward<T1>(contig_name),
        contig_begin_pos,
        std::forward<T2>(sequence_removed),
        std::forward<T3>(sequence_added),
        prior_model
    };
}
#endif /* defined(__Octopus__variant_factory__) */
