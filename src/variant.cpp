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
    return the_reference_allele_region_;
}

Allele Variant::get_reference_allele() const
{
    return Allele {the_reference_allele_region_, the_reference_allele_sequence_};
}

Allele Variant::get_alternative_allele() const
{
    return Allele {the_reference_allele_region_, the_alternative_allele_sequence_};
}

Variant::SizeType Variant::reference_allele_size() const noexcept
{
    return size(the_reference_allele_region_);
}

Variant::SizeType Variant::alternative_allele_size() const noexcept
{
    return static_cast<Variant::SizeType>(the_alternative_allele_sequence_.size());
}

const Variant::SequenceType& Variant::get_reference_allele_sequence() const noexcept
{
    return the_reference_allele_sequence_;
}

const Variant::SequenceType& Variant::get_alternative_allele_sequence() const noexcept
{
    return the_alternative_allele_sequence_;
}

std::ostream& operator<<(std::ostream& os, const Variant& a_variant)
{
    os << a_variant.get_region() << " " << a_variant.get_reference_allele_sequence() << " " <<
    a_variant.get_alternative_allele_sequence();
    return os;
}
