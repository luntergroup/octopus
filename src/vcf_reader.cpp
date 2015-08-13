//
//  vcf_reader.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_reader.h"

#include <fstream>

#include "vcf_header.h"
#include "vcf_record.h"

#include <iostream> // TEST

VcfReader::VcfReader(const fs::path& file_path)
:
file_path_ {file_path},
reader_ {file_path_, "r"}
{}

VcfHeader VcfReader::fetch_header()
{
    return reader_.fetch_header();
}

std::size_t VcfReader::num_records() const
{
    return reader_.num_records();
}

std::size_t VcfReader::num_records(const GenomicRegion& region) const
{
    return reader_.num_records(region);
}

std::vector<VcfRecord> VcfReader::fetch_records(Unpack level)
{
    return reader_.fetch_records((level == Unpack::All) ? HtslibBcfFacade::Unpack::All : HtslibBcfFacade::Unpack::AllButSamples);
}

std::vector<VcfRecord> VcfReader::fetch_records(const GenomicRegion& region, Unpack level)
{
    return reader_.fetch_records(region, (level == Unpack::All) ? HtslibBcfFacade::Unpack::All : HtslibBcfFacade::Unpack::AllButSamples);
}
