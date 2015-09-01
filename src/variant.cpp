//
//  variant.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variant.h"

const GenomicRegion& Variant::get_region() const noexcept
{
    return reference_allele_region_;
}

Allele Variant::get_reference_allele() const
{
    return Allele {reference_allele_region_, reference_allele_sequence_};
}

Allele Variant::get_alternative_allele() const
{
    return Allele {reference_allele_region_, alternative_allele_sequence_};
}

Variant::SizeType Variant::reference_allele_size() const noexcept
{
    return size(reference_allele_region_);
}

Variant::SizeType Variant::alternative_allele_size() const noexcept
{
    return static_cast<Variant::SizeType>(alternative_allele_sequence_.size());
}

const Variant::SequenceType& Variant::get_reference_allele_sequence() const noexcept
{
    return reference_allele_sequence_;
}

const Variant::SequenceType& Variant::get_alternative_allele_sequence() const noexcept
{
    return alternative_allele_sequence_;
}

bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_region() == rhs.get_region() &&
    lhs.get_reference_allele_sequence() == rhs.get_reference_allele_sequence() &&
    lhs.get_alternative_allele_sequence() == rhs.get_alternative_allele_sequence();
}

bool operator<(const Variant& lhs, const Variant& rhs)
{
    if (lhs.get_region() == rhs.get_region()) { // This check is required for consistency with operator==
        return (lhs.get_reference_allele_sequence() < rhs.get_reference_allele_sequence()) ? true :
        lhs.get_alternative_allele_sequence() < rhs.get_alternative_allele_sequence();
    } else {
        return lhs.get_region() < rhs.get_region();
    }
}

std::ostream& operator<<(std::ostream& os, const Variant& a_variant)
{
    os << a_variant.get_region() << " " << a_variant.get_reference_allele_sequence() << " " <<
        a_variant.get_alternative_allele_sequence();
    return os;
}
