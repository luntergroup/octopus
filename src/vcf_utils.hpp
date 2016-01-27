//
//  vcf_utils.hpp
//  Octopus
//
//  Created by Daniel Cooke on 11/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_utils__
#define __Octopus__vcf_utils__

#include <vector>
#include <string>
#include <map>

#include <boost/filesystem/path.hpp>

#include "vcf_type.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "vcf_reader.hpp"
#include "vcf_writer.hpp"

std::vector<std::string> get_contigs(const VcfHeader& header);

unsigned get_field_cardinality(const VcfHeader::KeyType& key, const VcfRecord& record);

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                                           const VcfHeader::KeyType& key);

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                             const VcfRecord::SampleIdType sample,
                                             const VcfHeader::KeyType& key);

void index_vcf(const boost::filesystem::path& vcf_path);
void index_vcf(const VcfReader& reader);
void index_vcfs(const std::vector<VcfReader>& readers);

std::vector<VcfReader> writers_to_readers(std::vector<VcfWriter>& writers);

VcfHeader merge(const std::vector<VcfHeader>& headers);

void merge(std::vector<VcfReader>& readers, const std::vector<std::string>& contigs, VcfWriter& result);

//VcfWriter merge(std::vector<VcfReader>& readers, boost::filesystem::path result_path);

#endif /* defined(__Octopus__vcf_utils__) */
