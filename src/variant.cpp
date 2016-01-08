//
//  variant.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant.hpp"

// public methods

Variant::Variant(const Allele& reference, const Allele& alternative)
:
reference_ {reference},
alternative_ {alternative}
{}

Variant::Variant(Allele&& reference, Allele&& alternative)
:
reference_ {std::move(reference)},
alternative_ {std::move(alternative)}
{}

const GenomicRegion& Variant::get_region() const noexcept
{
    return reference_.get_region();
}

const Allele& Variant::get_ref_allele() const noexcept
{
    return reference_;
}

const Allele& Variant::get_alt_allele() const noexcept
{
    return alternative_;
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
