//
//  vcf_utils.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "vcf_utils.hpp"

#include <algorithm>  // std::transform, std::adjacent_find
#include <functional> // std::not_equal_to
#include <stdexcept>

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
    
    VcfHeader::Builder hb {};
    
    hb.set_file_format(headers.front().get_file_format());
    hb.set_samples(headers.front().get_samples());
    
//    for (const auto& header : headers) {
//        
//    }
    
    return hb.build_once();
}

std::vector<VcfHeader> get_headers(const std::vector<VcfReader>& readers)
{
    std::vector<VcfHeader> result {};
    result.reserve(readers.size());
    
    std::transform(std::cbegin(readers), std::cend(readers), std::back_inserter(result),
                   [] (const auto& reader) { return reader.fetch_header(); });
    
    return result;
}

VcfWriter merge(const std::vector<VcfReader>& readers, fs::path result_path)
{
    auto header = merge(get_headers(readers));
    
    VcfWriter result {std::move(result_path), header};
    
    
    
    return result;
}
