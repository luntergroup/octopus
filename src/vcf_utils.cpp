//
//  vcf_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_utils.hpp"

#include <unordered_map>
#include <deque>
#include <queue>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <numeric>

#include "contig_region.hpp"
#include "genomic_region.hpp"

#include "htslib/vcf.h"
#include "htslib/tbx.h"

#include <iostream> // DEBUG

std::vector<std::string> get_contigs(const VcfHeader& header)
{
    std::vector<std::string> result {};
    
    const auto& contigs_fields = header.structured_fields("contig");
    
    result.reserve(contigs_fields.size());
    
    std::transform(std::cbegin(contigs_fields), std::cend(contigs_fields), std::back_inserter(result),
                   [] (const auto& field) { return field.at("ID"); });
    
    return result;
}

unsigned get_field_cardinality(const VcfHeader::StructuredKey& key, const VcfRecord& record)
{
    return 0;
}

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                                           const VcfHeader::StructuredKey& key)
{
    return get_typed_info_values(header, key, record.info_value(key.value));
}

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                             const VcfRecord::SampleIdType sample,
                                             const VcfHeader::StructuredKey& key)
{
    return get_typed_format_values(header, key, record.get_sample_value(sample, key.value));
}

void index_vcf(const boost::filesystem::path& vcf_path, const int lidx_shift)
{
    auto* const fp = hts_open(vcf_path.c_str(), "r");
    
    if (fp == nullptr) {
        throw std::ios::failure {vcf_path.c_str()};
    }
    
    const auto type = *hts_get_format(fp);
    
    hts_close(fp);
    
    if (type.format == bcf) {
        bcf_index_build(vcf_path.c_str(), lidx_shift);
    } else {
        tbx_index_build(vcf_path.c_str(), lidx_shift, &tbx_conf_vcf);
    }
}

void index_vcf(const VcfReader& reader, const int lidx_shift)
{
    index_vcf(reader.path(), lidx_shift);
}

void index_vcfs(const std::vector<VcfReader>& readers, const int lidx_shift)
{
    for (const auto& reader : readers) index_vcf(reader, lidx_shift);
}

std::vector<VcfReader> writers_to_readers(std::vector<VcfWriter>& writers)
{
    std::vector<VcfReader> result {};
    result.reserve(writers.size());
    
    for (auto& writer : writers) {
        auto path = writer.path();
        writer.close();
        result.emplace_back(std::move(path));
    }
    
    writers.clear();
    
    return result;
}

void copy(VcfReader& src, VcfWriter& dst)
{
    if (!dst.is_header_written()) {
        dst.write(src.fetch_header());
    }
    
    const auto records = src.fetch_records();
    
    write(records, dst);
}

bool all_same_format(const std::vector<VcfHeader>& headers)
{
    return std::adjacent_find(std::cbegin(headers), std::cend(headers),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.file_format() != rhs.file_format();
                              }) == std::cend(headers);
}

bool contain_same_samples(const std::vector<VcfHeader>& headers)
{
    return std::adjacent_find(std::cbegin(headers), std::cend(headers),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.samples() != rhs.samples();
                              }) == std::cend(headers);
}

bool all_equal(const std::vector<VcfHeader>& headers)
{
    return std::adjacent_find(std::cbegin(headers), std::cend(headers),
                              std::not_equal_to<VcfHeader>()) == std::cend(headers);
}

VcfHeader merge(const std::vector<VcfHeader>& headers)
{
    if (headers.empty()) return VcfHeader {};
    
    if (!all_same_format(headers)) {
        throw std::logic_error {"cannot merge VcfHeader's with different formats"};
    }
    
    if (!contain_same_samples(headers)) {
        throw std::logic_error {"cannot merge VcfHeader's that do not contain the same samples"};
    }
    
    if (all_equal(headers)) return headers.front();
    
    return headers.front(); // TODO: implement this
//    VcfHeader::Builder hb {};
//    
//    hb.set_file_format(headers.front().file_format());
//    hb.set_samples(headers.front().samples());
//    
////    for (const auto& header : headers) {
////        
////    }
//    
//    return hb.build_once();
}

std::vector<VcfHeader> get_headers(const std::vector<VcfReader>& readers)
{
    std::vector<VcfHeader> result {};
    result.reserve(readers.size());
    
    std::transform(std::cbegin(readers), std::cend(readers), std::back_inserter(result),
                   [] (const auto& reader) { return reader.fetch_header(); });
    
    return result;
}

using ReaderContigRecordCountMap = std::unordered_map<std::reference_wrapper<const VcfReader>,
                                                        std::unordered_map<std::string, size_t>>;

ReaderContigRecordCountMap get_contig_count_map(std::vector<VcfReader>& readers,
                                                const std::vector<std::string>& contigs)
{
    ReaderContigRecordCountMap result {};
    result.reserve(readers.size());
    
    for (auto& reader : readers) {
        std::unordered_map<std::string, size_t> contig_counts {};
        contig_counts.reserve(contigs.size());
        
        for (const auto& contig : contigs) {
            contig_counts.emplace(contig, reader.count_records(contig));
        }
        
        result.emplace(reader, std::move(contig_counts));
    }
    
    return result;
}

std::size_t calculate_total_number_of_records(const ReaderContigRecordCountMap& counts)
{
    return std::accumulate(std::cbegin(counts), std::cend(counts), std::size_t {0},
                           [] (const auto curr, const auto& reader_counts) {
                               const auto reader_total = std::accumulate(std::cbegin(reader_counts.second),
                                                                         std::cend(reader_counts.second),
                                                                         std::size_t {0},
                                                                         [] (const auto r_curr, const auto& p) {
                                                                             return r_curr + p.second;
                                                                         });
                               return curr + reader_total;
                           });
}

void merge(std::vector<VcfReader>& sources, VcfWriter& dst, const std::vector<std::string>& contigs)
{
    if (sources.empty()) return;
    
    if (sources.size() == 1) {
        copy(sources.front(), dst);
        return;
    }
    
    if (!dst.is_header_written()) {
        dst << merge(get_headers(sources));
    }
    
    auto reader_contig_counts = get_contig_count_map(sources, contigs);
    
    std::priority_queue<VcfRecord, std::deque<VcfRecord>, std::greater<VcfRecord>> record_queue {};
    
    auto total_record_count = calculate_total_number_of_records(reader_contig_counts);
    
    if (total_record_count <= 100000) {
        for (const auto& contig : contigs) {
            for (auto& reader : sources) {
                if (reader_contig_counts[reader].count(contig) == 1) {
                    auto records = reader.fetch_records(contig, VcfReader::UnpackPolicy::All);
                    for (auto&& record : records) {
                        record_queue.emplace(std::move(record));
                    }
                }
            }
            while (!record_queue.empty()) {
                dst << record_queue.top();
                record_queue.pop();
            }
        }
        return;
    }
    
    static constexpr GenomicRegion::SizeType buffer_size {10000}; // maybe make contig dependent
    
    for (const auto& contig : contigs) {
        GenomicRegion region {contig, 0, buffer_size};
        
        bool all_done {false};
        
        while (!all_done) {
            all_done = true;
            
            for (auto& reader : sources) {
                if (reader_contig_counts[reader].count(contig) == 1) {
                    auto records = reader.fetch_records(region, VcfReader::UnpackPolicy::All);
                    
                    reader_contig_counts[reader][contig] -= records.size();
                    
                    if (reader_contig_counts[reader][contig] > 0) all_done = false;
                    
                    for (auto&& record : records) {
                        record_queue.emplace(std::move(record));
                    }
                }
            }
            
            while (!record_queue.empty()) {
                dst << record_queue.top();
                record_queue.pop();
            }
            
            region = shift(region, buffer_size);
        }
    }
}

void merge(std::vector<VcfReader>& sources, VcfWriter& dst)
{
    if (sources.empty()) return;
    
    if (sources.size() == 1) {
        copy(sources.front(), dst);
        return;
    }
    
    const auto header = merge(get_headers(sources));
    
    if (!dst.is_header_written()) {
        dst << header;
    }
    
    const auto contigs = get_contigs(header);
    
    return merge(sources, dst, contigs);
}
