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

#include <boost/filesystem/path.hpp>

#include "vcf_type.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "vcf_reader.hpp"
#include "vcf_writer.hpp"

std::vector<std::string> get_contigs(const VcfHeader& header);

unsigned get_field_cardinality(const VcfHeader::StructuredKey& key, const VcfRecord& record);

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                                           const VcfHeader::StructuredKey& key);

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                             const VcfRecord::SampleIdType sample,
                                             const VcfHeader::StructuredKey& key);

void index_vcf(const boost::filesystem::path& vcf_path, int lidx_shift = 14);
void index_vcf(const VcfReader& reader, int lidx_shift = 14);
void index_vcfs(const std::vector<VcfReader>& readers, int lidx_shift = 14);

std::vector<VcfReader> writers_to_readers(std::vector<VcfWriter>& writers);

void copy(const VcfReader& src, VcfWriter& dst);

VcfHeader merge(const std::vector<VcfHeader>& headers);

void merge(const std::vector<VcfReader>& sources, VcfWriter& dst, const std::vector<std::string>& contigs);
void merge(const std::vector<VcfReader>& sources, VcfWriter& dst);

#endif /* defined(__Octopus__vcf_utils__) */
