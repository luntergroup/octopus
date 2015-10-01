//
//  vcf_reader.cpp
//  Octopus
//
//  Created by Daniel Cooke on 28/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_reader.hpp"

#include <fstream>
#include <iostream>
#include <stdexcept>

#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "htslib_bcf_facade.hpp"
#include "vcf_parser.hpp"

std::unique_ptr<IVcfReaderImpl> make_vcf_reader(const fs::path& file_path)
{
    auto file_type = file_path.extension().string();
    
    if (file_type == ".vcf") {
        auto vcf_file_size = fs::file_size(file_path);
        if (vcf_file_size > 1e9) { // 1GB
            throw std::runtime_error {"VCF file " + file_path.string() + " is too big"};
        }
        
        return std::make_unique<VcfParser>(file_path);
    } else {
        return std::make_unique<HtslibBcfFacade>(file_path, "r");
    }
}

VcfReader::VcfReader(const fs::path& file_path)
:
file_path_ {file_path},
reader_ {make_vcf_reader(file_path)}
{}

const fs::path VcfReader::path() const
{
    return file_path_;
}

VcfHeader VcfReader::fetch_header()
{
    return reader_->fetch_header();
}

size_t VcfReader::count_records() const
{
    return reader_->count_records();
}

size_t VcfReader::count_records(const GenomicRegion& region) const
{
    return reader_->count_records(region);
}

std::vector<VcfRecord> VcfReader::fetch_records(Unpack level)
{
    return reader_->fetch_records((level == Unpack::All) ? HtslibBcfFacade::Unpack::All : HtslibBcfFacade::Unpack::AllButSamples);
}

std::vector<VcfRecord> VcfReader::fetch_records(const GenomicRegion& region, Unpack level)
{
    return reader_->fetch_records(region, (level == Unpack::All) ? HtslibBcfFacade::Unpack::All : HtslibBcfFacade::Unpack::AllButSamples);
}
