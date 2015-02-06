//
//  variant.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_h
#define Octopus_variant_h

#include <string>
#include <memory>
#include <cstdint>
#include <functional>

#include "genome_region.h"
#include "variant_prior_model.h"

using std::size_t;

/*
 Class representing all variant types.
 
 The
 
 */
class Variant
{
public:
    Variant() = delete;
    Variant(GenomeRegion ref_region, std::string sequence_added, std::string sequence_removed,
            VariantPriorModel* prior_model);
    ~Variant();
    
    GenomeRegion get_ref_region() const noexcept;
    unsigned long get_num_supporting_reads() const noexcept;
    bool overlaps(const Variant& other) const noexcept;
    
    double get_prior_probability() const noexcept;
    
private:
    GenomeRegion ref_region_;
    std::string sequence_added_, sequence_removed_;
    
    unsigned long num_supporting_reads_;
    std::shared_ptr<VariantPriorModel> prior_model_;
};

inline bool operator==(const Variant& lhs, const Variant& rhs);
inline bool operator!=(const Variant& lhs, const Variant& rhs);
inline bool operator< (const Variant& lhs, const Variant& rhs);
inline bool operator> (const Variant& lhs, const Variant& rhs);
inline bool operator<=(const Variant& lhs, const Variant& rhs);
inline bool operator>=(const Variant& lhs, const Variant& rhs);

namespace std {
    template <> struct hash<Variant>
    {
        size_t operator()(const Variant& v) const
        {
            return hash<size_t>()(v.get_ref_region().begin); //TODO: do something better!
        }
    };
}

#endif
