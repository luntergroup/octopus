//
//  variant.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant.hpp"

Variant::Variant(Allele reference, Allele alternative)
:
ref_ {std::move(reference)},
alt_ {std::move(alternative)}
{}

const Allele& Variant::get_ref_allele() const noexcept
{
    return ref_;
}

const Allele&  Variant::get_alt_allele() const noexcept
{
    return alt_;
}

// non-member methods

const Variant::SequenceType& get_ref_sequence(const Variant& variant)
{
    return variant.get_ref_allele().get_sequence();
}

const Variant::SequenceType& get_alt_sequence(const Variant& variant)
{
    return variant.get_alt_allele().get_sequence();
}

Variant::SizeType ref_sequence_size(const Variant& variant)
{
    return sequence_size(variant.get_ref_allele());
}

Variant::SizeType alt_sequence_size(const Variant& variant)
{
    return sequence_size(variant.get_alt_allele());
}

bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_ref_allele() == rhs.get_ref_allele() && get_alt_sequence(lhs) == get_alt_sequence(rhs);
}

bool operator<(const Variant& lhs, const Variant& rhs)
{
    return (lhs.get_ref_allele() < rhs.get_ref_allele() ||
            (lhs.get_ref_allele() == rhs.get_ref_allele() && lhs.get_alt_allele() < rhs.get_alt_allele()));
}

std::ostream& operator<<(std::ostream& os, const Variant& variant)
{
    os << get_region(variant) << " " << get_ref_sequence(variant) << " " << get_alt_sequence(variant);
    return os;
}
