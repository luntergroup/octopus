// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "reference_call.hpp"

#include "concepts/mappable.hpp"

namespace octopus {

const GenomicRegion& ReferenceCall::mapped_region() const noexcept
{
    return octopus::mapped_region(reference_);
}

const Allele& ReferenceCall::reference() const noexcept
{
    return reference_;
}

bool ReferenceCall::is_represented(const Allele& allele) const noexcept
{
    return allele == reference_;
}

void ReferenceCall::replace_called_alleles(const char old_base, const char replacement_base)
{
    const auto& ref_sequence = reference_.sequence();
    
    auto it = std::find(std::cbegin(ref_sequence), std::cend(ref_sequence), old_base);
    if (it != std::cend(ref_sequence)) {
        Allele::NucleotideSequence new_sequence {};
        new_sequence.reserve(ref_sequence.size());
        new_sequence.insert(std::end(new_sequence), std::cbegin(ref_sequence), it);
        std::replace_copy(it, std::cend(ref_sequence), std::back_inserter(new_sequence),
                          old_base, replacement_base);
        reference_ = Allele {reference_.mapped_region(), std::move(new_sequence)};
    }
}

void ReferenceCall::replace(const Allele& old, Allele replacement)
{
    if (reference_ == old) reference_ = std::move(replacement);
}

void ReferenceCall::replace_uncalled_genotype_alleles(const Allele& replacement, const char ignore)
{
    
}

void ReferenceCall::decorate(VcfRecord::Builder& record) const
{
    
}

std::unique_ptr<Call> ReferenceCall::do_clone() const
{
    return std::make_unique<ReferenceCall>(*this);
}

} // namespace octopus
