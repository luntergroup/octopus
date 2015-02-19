//
//  htslib_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_facade.h"

#include <sstream>

HtslibFacade::HtslibFacade(const std::string& htslib_file_path)
:the_filename_ {htslib_file_path},
 the_file_ {hts_open(the_filename_.c_str(), "r"), htslib_file_deleter},
 the_header_ {sam_hdr_read(the_file_.get()), htslib_header_deleter},
 the_index_ {sam_index_load(the_file_.get(), the_filename_.c_str()), htslib_index_deleter},
 contig_name_map_ {},
 sample_id_map_ {}
{
    contig_name_map_ = get_htslib_tid_to_contig_name_map();
    sample_id_map_   = get_read_group_to_sample_id_map();
}

void HtslibFacade::close()
{
    // TODO: what should this do?
}

std::vector<std::string> HtslibFacade::get_sample_ids()
{
    std::vector<std::string> result {};
    for (const auto pair : sample_id_map_) {
        if (std::find(std::cbegin(result), std::cend(result), pair.second) == std::cend(result)) {
            result.emplace_back(pair.second);
        }
    }
    return result;
}

std::vector<std::string> HtslibFacade::get_read_groups_in_sample(const std::string& a_sample_id)
{
    std::vector<std::string> result {};
    for (const auto pair : sample_id_map_) {
        if (pair.second == a_sample_id) result.emplace_back(pair.first);
    }
    return result;
}

HtslibFacade::SampleIdToReadsMap HtslibFacade::fetch_reads(const GenomicRegion& a_region)
{
    HtslibIterator it {*this, a_region};
    SampleIdToReadsMap the_reads {};
    while (++it) {
        try {
            auto&& a_read_and_its_group = *it;
            auto the_sample_id = sample_id_map_[a_read_and_its_group.second];
            the_reads[std::move(the_sample_id)].emplace_back(std::move(a_read_and_its_group.first));
        } catch (const std::runtime_error& e) {
            // TODO: log maybe?
            // There isn't much we can do here
        }
    }
    return the_reads;
}

std::vector<std::string> HtslibFacade::get_reference_contig_names()
{
    std::vector<std::string> result {};
    result.reserve(get_num_reference_contigs());
    for (uint_fast32_t i {0}; i < get_num_reference_contigs(); ++i) {
        result.emplace_back(the_header_->target_name[i]);
    }
    return result;
}

std::vector<GenomicRegion> HtslibFacade::get_regions_in_file()
{
    std::vector<GenomicRegion> result {};
    for (uint_fast32_t i {0}; i < get_num_reference_contigs(); ++i) {
        auto contig_name = get_reference_contig_name(i);
        // CRAM files don't seem to have the same index stats as BAM files so
        // we don't know which contigs have been mapped to
        if (the_file_->is_cram || get_num_mapped_reads(contig_name) > 0) {
            result.emplace_back(contig_name, 0, get_reference_contig_size(contig_name));
        }
    }
    return result;
}

HtslibFacade::ReadGroupToSampleIdMap HtslibFacade::get_read_group_to_sample_id_map() const
{
    std::string the_header_text (the_header_->text, the_header_->l_text);
    std::stringstream ss {the_header_text};
    std::string line {};
    ReadGroupToSampleIdMap result {};
    unsigned num_read_groups {0};
    while (std::getline(ss, line, '\n')) {
        if (is_type(line, Read_group_tag)) {
            if (!has_tag(line, Read_group_id_tag) || !has_tag(line, Sample_id_tag)) {
                // The SAM specification does not mark the sample id tag 'SM' as a required
                // field, however we can't do much without it.
                throw std::runtime_error {"bad read file"};
            }
            result.emplace(std::make_pair(get_tag_value(line, Read_group_id_tag),
                                                        get_tag_value(line, Sample_id_tag)));
            ++num_read_groups;
        }
    }
    if (num_read_groups == 0) throw std::runtime_error {"bad read file"};
    return result;
}

HtslibFacade::HtsTidToContigNameMap HtslibFacade::get_htslib_tid_to_contig_name_map() const
{
    HtsTidToContigNameMap result {};
    result.reserve(the_header_->n_targets);
    for (uint_fast32_t i {0}; i < the_header_->n_targets; ++i) {
        result.emplace(std::make_pair(i, the_header_->target_name[i]));
    }
    return result;
}

HtslibFacade::HtslibIterator::HtslibIterator(HtslibFacade& hts_facade, const GenomicRegion& a_region)
:hts_facade_ {hts_facade},
 the_iterator_ {sam_itr_querys(hts_facade_.the_index_.get(), hts_facade_.the_header_.get(),
                               to_string(a_region).c_str()), htslib_iterator_deleter},
 the_bam1_ {bam_init1(), htslib_bam1_deleter}
{}

std::pair<AlignedRead, std::string> HtslibFacade::HtslibIterator::operator*() const
{
    auto the_qualities = get_qualities();
    if (the_qualities.empty() || the_qualities[0] == 0xff) throw std::runtime_error {"bad sequence"};
    auto the_cigar_string = get_cigar_string();
    if (the_cigar_string.empty()) throw std::runtime_error {"bad sequence"};
    auto c = the_bam1_->core;
    auto read_start = static_cast<uint_fast32_t>(get_soft_clipped_read_begin(the_cigar_string, c.pos));
    
    return {AlignedRead {
        GenomicRegion(get_contig_name(c.tid), read_start, read_start + get_sequence_length()),
        get_sequence(),
        std::move(the_qualities),
        std::move(the_cigar_string),
        static_cast<uint_fast32_t>(c.isize),
        get_contig_name(c.mtid),
        static_cast<uint_fast32_t>(c.mpos),
        static_cast<uint_fast8_t>(c.qual)
    }, get_read_group()};
}

int32_t HtslibFacade::get_htslib_tid(const std::string& reference_contig_name) const
{
    // Could avoid this lookup if stored inverse tid mappings, but I don't think this
    // function will be called much, so probably better not to explicitly store map.
    for (const auto& pair : contig_name_map_) {
        if (pair.second == reference_contig_name) return pair.first;
    }
    throw std::runtime_error {"reference contig not found"};
}
