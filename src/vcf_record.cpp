//
//  vcf_record.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_record.h"

#include <algorithm> // std::find
#include <iterator>  // std::cbegin, std::cend

const std::string& VcfRecord::get_chrom() const noexcept
{
    return chrom_;
}

VcfRecord::SizeType VcfRecord::get_pos() const noexcept
{
    return pos_;
}

const std::string& VcfRecord::get_id() const noexcept
{
    return id_;
}

const VcfRecord::SequenceType& VcfRecord::get_ref() const noexcept
{
    return ref_;
}

unsigned VcfRecord::get_num_alt_alleles() const noexcept
{
    return static_cast<unsigned>(alt_.size());
}

const VcfRecord::SequenceType& VcfRecord::get_alt(unsigned n) const noexcept
{
    return alt_.at(n);
}

VcfRecord::QualityType VcfRecord::get_qual() const noexcept
{
    return qual_;
}

bool VcfRecord::has_filter(const std::string& filter) const noexcept
{
    return std::find(std::cbegin(filter_), std::cend(filter_), filter) != std::cend(filter_);
}

unsigned VcfRecord::get_num_samples() const noexcept
{
    return static_cast<unsigned>(format_.size());
}

// non-member functions

std::ostream& operator<<(std::ostream& os, const VcfRecord& record)
{
    os << record.get_chrom() << "\t" << record.get_pos() << "\t" << record.get_id() << "\t" << record.get_ref() << "\t";
    if (record.get_num_alt_alleles() == 0) {
        os << ".\t";
    } else {
        for (unsigned i {}; i < record.get_num_alt_alleles() - 1; ++i) {
            os << record.get_alt(i) << ",";
        }
        os << record.get_alt(record.get_num_alt_alleles() - 1) << "\t";
    }
    os << static_cast<unsigned>(record.get_qual());
    return os;
}
