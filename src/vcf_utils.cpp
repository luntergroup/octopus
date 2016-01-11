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

#include "contig_region.hpp"
#include "genomic_region.hpp"

#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "htslib/tbx.h"

#include <iostream> // DEBUG

std::vector<std::string> get_contigs(const VcfHeader& header)
{
    std::vector<std::string> result {};
    
    const auto& contigs_fields = header.get_structured_fields("contig");
    
    result.reserve(contigs_fields.size());
    
    std::transform(std::cbegin(contigs_fields), std::cend(contigs_fields), std::back_inserter(result),
                   [] (const auto& field) {
                       return field.at("ID");
                   });
    
    return result;
}

unsigned get_field_cardinality(const VcfHeader::KeyType& key, const VcfRecord& record)
{
    return 0;
}

std::vector<VcfType> get_typed_info_values(const VcfHeader& header, const VcfRecord& record,
                              const VcfHeader::KeyType& key)
{
    return get_typed_info_values(header, key, record.get_info_value(key));
}

std::vector<VcfType> get_typed_format_values(const VcfHeader& header, const VcfRecord& record,
                                const VcfRecord::SampleIdType sample, const VcfHeader::KeyType& key)
{
    return get_typed_format_values(header, key, record.get_sample_value(sample, key));
}

void index_vcf(const boost::filesystem::path& vcf_file)
{
    
}

void index_vcf(const boost::filesystem::path& vcf_file, const boost::filesystem::path& out_index_path)
{
    htsFile *fp {hts_open(vcf_file.c_str(), "r")};
    htsFormat type {*hts_get_format(fp)};
    hts_close(fp);
    
//    if ((type.format != .bcf && type.format ! =vcf) || type.compression != bgzf) {
//        
//    }
}

bool all_same_format(const std::vector<VcfHeader>& headers)
{
    using std::cbegin; using std::cend;
    return std::adjacent_find(cbegin(headers), cend(headers),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.get_file_format() != rhs.get_file_format();
                              }) == cend(headers);
}

bool contain_same_samples(const std::vector<VcfHeader>& headers)
{
    using std::cbegin; using std::cend;
    return std::adjacent_find(cbegin(headers), cend(headers),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.get_samples() != rhs.get_samples();
                              }) == cend(headers);
}

bool all_the_same(const std::vector<VcfHeader>& headers)
{
    using std::cbegin; using std::cend;
    return std::adjacent_find(cbegin(headers), cend(headers), std::not_equal_to<VcfHeader>()) == cend(headers);
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
    
    if (all_the_same(headers)) return headers.front();
    
    
    return headers.front(); // TODO: implement this
//    VcfHeader::Builder hb {};
//    
//    hb.set_file_format(headers.front().get_file_format());
//    hb.set_samples(headers.front().get_samples());
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

auto get_contig_count_map(std::vector<VcfReader>& readers, const std::vector<std::string>& contigs)
{
    std::unordered_map<std::reference_wrapper<const VcfReader>, std::unordered_map<std::string, size_t>> result {};
    
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

VcfWriter merge(std::vector<VcfReader>& readers, boost::filesystem::path result_path)
{
    auto header = merge(get_headers(readers));
    
    auto contigs = get_contigs(header);
    
    std::sort(std::begin(contigs), std::end(contigs));
    
    auto reader_contig_counts = get_contig_count_map(readers, contigs);
    
    std::priority_queue<VcfRecord, std::deque<VcfRecord>, std::greater<VcfRecord>> record_queue {};
    
    static constexpr GenomicRegion::SizeType buffer_size {10000}; // maybe make contig dependent
    
    VcfWriter result {std::move(result_path), header};
    
    for (const auto& contig : contigs) {
        GenomicRegion region {contig, 0, buffer_size};
        
        bool all_done {false};
        
        while (!all_done) {
            all_done = true;
            
            for (auto& reader : readers) {
                if (reader_contig_counts[reader].count(contig) == 1) {
                    auto records = reader.fetch_records(region, VcfReader::Unpack::All);
                    
                    reader_contig_counts[reader][contig] -= records.size();
                    
                    if (reader_contig_counts[reader][contig] > 0) all_done = false;
                    
                    for (auto&& record : records) {
                        record_queue.emplace(std::move(record));
                    }
                }
            }
            
            while (!record_queue.empty()) {
                result.write(record_queue.top());
                record_queue.pop();
            }
            
            region = shift(region, buffer_size);
        }
    }
    
    return result;
}
