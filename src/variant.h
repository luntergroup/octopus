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
#include <cstdint>
#include <functional>

#include "genome_region.h"

using std::size_t;

/*
 Class representing all variant types.
 
 The
 
 */
class Variant
{
public:
    Variant() = delete;
    Variant(GenomicRegion ref_region, std::string sequence_added, std::string sequence_removed,
            std::function<double()> prior_model);
    ~Variant();
    
    GenomicRegion get_ref_region() const noexcept;
    const std::string& get_sequence_added() const noexcept;
    const std::string& get_sequence_removed() const noexcept;
    unsigned long get_num_supporting_reads() const noexcept;
    
    void add_support(unsigned long num_reads) noexcept;
    bool overlaps(const Variant& other) const noexcept;
    
    double get_prior_probability() const noexcept;
    
private:
    GenomicRegion ref_region_;
    std::string sequence_added_, sequence_removed_;
    
    unsigned long num_supporting_reads_;
    std::function<double()> prior_model_;
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
