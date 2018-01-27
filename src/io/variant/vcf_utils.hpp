// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef vcf_utils_hpp
#define vcf_utils_hpp

#include <vector>
#include <string>

#include <boost/filesystem/path.hpp>

#include "vcf_type.hpp"
#include "vcf_header.hpp"
#include "vcf_record.hpp"
#include "vcf_reader.hpp"
#include "vcf_writer.hpp"

namespace octopus {

std::vector<std::string> get_contigs(const VcfHeader& header);

unsigned get_field_cardinality(const VcfHeader::StructuredKey& key, const VcfRecord& record);

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                                           const VcfHeader::StructuredKey& key);

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                             const VcfRecord::SampleName sample,
                                             const VcfHeader::StructuredKey& key);

bool is_indexable(const boost::filesystem::path& vcf_path);

void index_vcf(const boost::filesystem::path& vcf_path);
void index_vcf(const VcfReader& reader);
void index_vcfs(const std::vector<VcfReader>& readers);

std::vector<VcfReader> writers_to_readers(std::vector<VcfWriter>&& writers);

void copy(const VcfReader& src, VcfWriter& dst);

void sort(const VcfReader& src, VcfWriter& dst);

VcfHeader merge(const std::vector<VcfHeader>& headers);

void merge(const std::vector<VcfReader>& sources, VcfWriter& dst,
           const std::vector<std::string>& contigs);

void merge(const std::vector<VcfReader>& sources, VcfWriter& dst);

void convert_to_legacy(const VcfReader& src, VcfWriter& dst, bool remove_ref_pad_duplicates = true);

} // namespace octopus    

#endif
