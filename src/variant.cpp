//
//  variant.cpp
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant.h"

Variant::Variant(GenomicRegion ref_region, std::string sequence_added, std::string sequence_removed,
                 std::function<double()> prior_model)
: ref_region_(ref_region), sequence_added_(sequence_added), sequence_removed_(sequence_removed),
    prior_model_(prior_model)
{
    // TODO: Other setup chores
}

Variant::~Variant() {} // Nothing allocated

GenomicRegion Variant::get_ref_region() const noexcept
{
    return ref_region_;
}

const std::string& Variant::get_sequence_added() const noexcept
{
    return sequence_added_;
}

const std::string& Variant::get_sequence_removed() const noexcept
{
    return sequence_removed_;
}

unsigned long Variant::get_num_supporting_reads() const noexcept
{
    return num_supporting_reads_;
}

void Variant::add_support(unsigned long num_reads) noexcept
{
    num_supporting_reads_ += num_reads;
}

bool Variant::overlaps(const Variant& other) const noexcept
{
    return false; // TODO: do actual comparison
}

double Variant::get_prior_probability() const noexcept
{
    return prior_model_();
}

inline bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_ref_region() == rhs.get_ref_region() &&
        lhs.get_sequence_added() == rhs.get_sequence_added() &&
        lhs.get_sequence_removed() == rhs.get_sequence_removed();
}
inline bool operator!=(const Variant& lhs, const Variant& rhs) {return !operator==(lhs,rhs);}
inline bool operator< (const Variant& lhs, const Variant& rhs)
{
    return lhs.get_ref_region() < rhs.get_ref_region();
}
inline bool operator> (const Variant& lhs, const Variant& rhs) {return  operator< (rhs,lhs);}
inline bool operator<=(const Variant& lhs, const Variant& rhs) {return !operator> (lhs,rhs);}
inline bool operator>=(const Variant& lhs, const Variant& rhs) {return !operator< (lhs,rhs);}
