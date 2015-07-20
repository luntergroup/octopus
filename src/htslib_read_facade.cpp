//
//  htslib_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_read_facade.h"

#include <sstream>
#include <cmath>   // std::abs
#include <utility> // std::move

#include "cigar_string.h"

#include <iostream> // TEST

HtslibReadFacade::HtslibReadFacade(const fs::path& file_path)
:
file_path_ {file_path},
hts_file_ {hts_open(file_path_.string().c_str(), "r"), htslib_file_deleter},
hts_header_ {sam_hdr_read(hts_file_.get()), htslib_header_deleter},
hts_index_ {sam_index_load(hts_file_.get(), file_path_.string().c_str()), htslib_index_deleter},
hts_tid_map_ {},
contig_name_map_ {},
sample_id_map_ {}
{
    if (hts_file_ == nullptr) {
        throw std::runtime_error {"could not open " + file_path_.string()};
    }
    
    if (hts_header_ == nullptr) {
        throw std::runtime_error {"could not open file header for " + file_path_.string()};
    }
    
    if (hts_index_ == nullptr) {
        throw std::runtime_error {"could not open index file for " + file_path_.string()};
    }
    
    std::tie(hts_tid_map_, contig_name_map_) = get_htslib_tid_maps();
    sample_id_map_ = get_read_group_to_sample_id_map();
}

void HtslibReadFacade::close()
{
    // TODO: what should this do?
}

unsigned HtslibReadFacade::get_num_reference_contigs() noexcept
{
    return hts_header_->n_targets;
}

HtslibReadFacade::SizeType HtslibReadFacade::get_reference_contig_size(const std::string& contig_name)
{
    return hts_header_->target_len[get_htslib_tid(contig_name)];
}

uint64_t HtslibReadFacade::get_num_mapped_reads(const std::string& reference_contig_name) const
{
    uint64_t num_mapped, num_unmapped;
    hts_idx_get_stat(hts_index_.get(), get_htslib_tid(reference_contig_name), &num_mapped, &num_unmapped);
    return num_mapped;
}

std::vector<HtslibReadFacade::SampleIdType> HtslibReadFacade::get_sample_ids()
{
    std::vector<HtslibReadFacade::SampleIdType> result {};
    for (const auto pair : sample_id_map_) {
        if (std::find(std::cbegin(result), std::cend(result), pair.second) == std::cend(result)) {
            result.emplace_back(pair.second);
        }
    }
    return result;
}

std::vector<std::string> HtslibReadFacade::get_read_groups_in_sample(const std::string& a_sample_id)
{
    std::vector<std::string> result {};
    for (const auto pair : sample_id_map_) {
        if (pair.second == a_sample_id) result.emplace_back(pair.first);
    }
    return result;
}

std::size_t HtslibReadFacade::get_num_reads(const GenomicRegion& a_region)
{
    std::size_t result {0};
    HtslibIterator it {*this, a_region};
    while (++it) { ++result; }
    return result;
}

HtslibReadFacade::SampleIdToReadsMap HtslibReadFacade::fetch_reads(const GenomicRegion& a_region)
{
    HtslibIterator it {*this, a_region};
    SampleIdToReadsMap the_reads {};
    
    while (++it) {
        try {
            auto a_read_and_its_group = *it;
            const auto& the_sample_id = sample_id_map_.at(a_read_and_its_group.second);
            the_reads[std::move(the_sample_id)].emplace_back(std::move(a_read_and_its_group.first));
        } catch (const std::runtime_error& e) {
            // TODO: log maybe?
            // There isn't much we can do here
            //std::cout << e.what() << std::endl;
        }
    }
    
    for (auto& sample_reads_pair : the_reads) {
        sample_reads_pair.second.shrink_to_fit();
    }
    
    return the_reads;
}

std::vector<std::string> HtslibReadFacade::get_reference_contig_names()
{
    std::vector<std::string> result {};
    result.reserve(get_num_reference_contigs());
    
    for (HtsTidType hts_tid {0}; hts_tid < get_num_reference_contigs(); ++hts_tid) {
        result.emplace_back(get_reference_contig_name(hts_tid));
    }
    
    return result;
}

std::vector<GenomicRegion> HtslibReadFacade::get_possible_regions_in_file()
{
    std::vector<GenomicRegion> result {};
    result.reserve(get_num_reference_contigs());
    std::string contig_name;
    
    for (HtsTidType hts_tid {0}; hts_tid < get_num_reference_contigs(); ++hts_tid) {
        contig_name = get_reference_contig_name(hts_tid);
        // CRAM files don't seem to have the same index stats as BAM files so
        // we don't know which contigs have been mapped to
        if (hts_file_->is_cram || get_num_mapped_reads(contig_name) > 0) {
            result.emplace_back(contig_name, 0, get_reference_contig_size(contig_name));
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

HtslibReadFacade::ReadGroupToSampleIdMap HtslibReadFacade::get_read_group_to_sample_id_map() const
{
    std::string the_header_text (hts_header_->text, hts_header_->l_text);
    std::stringstream ss {the_header_text};
    std::string line {};
    ReadGroupToSampleIdMap result {};
    unsigned num_read_groups {0};
    
    while (std::getline(ss, line, '\n')) {
        if (is_tag_type(line, Read_group_tag)) {
            if (!(has_tag(line, Read_group_id_tag) && has_tag(line, Sample_id_tag))) {
                // The SAM specification does not specify the sample id tag 'SM' as a required
                // field, however we can't do much without it.
                throw std::runtime_error {"no read group tag in sample (SM) header line in file " +
                                           file_path_.string()};
            }
            result.emplace(get_tag_value(line, Read_group_id_tag), get_tag_value(line, Sample_id_tag));
            ++num_read_groups;
        }
    }
    
    if (num_read_groups == 0) {
        throw std::runtime_error {"no read groups found in file " + file_path_.string()};
    }
    
    return result;
}

std::pair<HtslibReadFacade::ContigNameToHtsTidMap, HtslibReadFacade::HtsTidToContigNameMap>
HtslibReadFacade::get_htslib_tid_maps()
{
    ContigNameToHtsTidMap hts_tid_map_ {};
    HtsTidToContigNameMap result {};
    
    hts_tid_map_.reserve(get_num_reference_contigs());
    result.reserve(get_num_reference_contigs());
    
    for (HtsTidType hts_tid {0}; hts_tid < get_num_reference_contigs(); ++hts_tid) {
        hts_tid_map_.emplace(hts_header_->target_name[hts_tid], hts_tid);
        result.emplace(hts_tid, hts_header_->target_name[hts_tid]);
    }
    
    return {hts_tid_map_, result};
}

std::string HtslibReadFacade::get_reference_contig_name(HtsTidType hts_tid) const
{
    return std::string(hts_header_->target_name[hts_tid]);
}

bool HtslibReadFacade::is_tag_type(const std::string& header_line, const char* tag) const
{
    return header_line.compare(1, 2, tag) == 0;
}

bool HtslibReadFacade::has_tag(const std::string& header_line, const char* tag) const
{
    return header_line.find(tag) != std::string::npos;
}

std::string HtslibReadFacade::get_tag_value(const std::string& line, const char* tag) const
{
    // format TAG:VALUE\t
    auto tag_position = line.find(tag);
    if (tag_position != std::string::npos) {
        auto value_position = line.find(':', tag_position) + 1;
        auto tag_value_size = line.find('\t', value_position) - value_position;
        return line.substr(value_position, tag_value_size);
    }
    
    throw std::runtime_error {"no " + std::string {tag} + " tag"};
}

HtslibReadFacade::HtsTidType HtslibReadFacade::get_htslib_tid(const std::string& reference_contig_name) const
{
    if (hts_tid_map_.count(reference_contig_name) == 0) {
        throw std::runtime_error {"reference contig not found"};
    }
    
    return hts_tid_map_.at(reference_contig_name);
}

//
// HtslibIterator
//

HtslibReadFacade::HtslibIterator::HtslibIterator(HtslibReadFacade& hts_facade, const GenomicRegion& a_region)
:
hts_facade_ {hts_facade},
hts_iterator_ {sam_itr_querys(hts_facade_.hts_index_.get(), hts_facade_.hts_header_.get(),
                              to_string(a_region).c_str()), htslib_iterator_deleter},
hts_bam1_ {bam_init1(), htslib_bam1_deleter}
{
    if (hts_iterator_ == nullptr) {
        throw std::runtime_error {"could not load read iterator for " + hts_facade.file_path_.string()};
    }
    
    if (hts_bam1_ == nullptr) {
        throw std::runtime_error {"error creating bam1 for " + hts_facade.file_path_.string()};
    }
}

bool HtslibReadFacade::HtslibIterator::operator++()
{
    return sam_itr_next(hts_facade_.hts_file_.get(), hts_iterator_.get(), hts_bam1_.get()) >= 0;
}

std::pair<AlignedRead, HtslibReadFacade::SampleIdType> HtslibReadFacade::HtslibIterator::operator*() const
{
    auto the_qualities = get_qualities();
    
    if (the_qualities.empty() || the_qualities[0] == 0xff) {
        throw std::runtime_error {"improper sequence data in read " + get_read_name() +
                                    " in file " + hts_facade_.file_path_.string()};
    }
    
    auto the_cigar_string = make_cigar_string();
    
    if (the_cigar_string.empty()) {
        throw std::runtime_error {"improper cigar data in read " + get_read_name() +
                                    " in file " + hts_facade_.file_path_.string()};
    }
    
    auto c = hts_bam1_->core;
    
    auto read_start = static_cast<AlignedRead::SizeType>(soft_clipped_read_begin(the_cigar_string, c.pos));
    
    if (c.mtid == -1) { // TODO: check if this is always true
        return {AlignedRead {
            GenomicRegion {get_contig_name(c.tid), read_start, read_start +
                reference_size<GenomicRegion::SizeType>(the_cigar_string)},
            get_sequence(),
            std::move(the_qualities),
            std::move(the_cigar_string),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags()
        }, get_read_group()};
    } else {
        return {AlignedRead {
            GenomicRegion {get_contig_name(c.tid), read_start, read_start +
                reference_size<GenomicRegion::SizeType>(the_cigar_string)},
            get_sequence(),
            std::move(the_qualities),
            std::move(the_cigar_string),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags(),
            get_contig_name(c.mtid),
            static_cast<AlignedRead::SizeType>(c.mpos),
            static_cast<AlignedRead::SizeType>(std::abs(c.isize)),
            get_next_segment_flags()
        }, get_read_group()};
    }
}

HtslibReadFacade::SizeType HtslibReadFacade::HtslibIterator::get_read_start() const noexcept
{
    return static_cast<HtslibReadFacade::SizeType>(hts_bam1_->core.pos);
}

uint32_t HtslibReadFacade::HtslibIterator::get_sequence_length() const noexcept
{
    return hts_bam1_->core.l_qseq;
}

char HtslibReadFacade::HtslibIterator::get_base(uint8_t* a_htslib_sequence, uint32_t index) const noexcept
{
    static constexpr const char* symbol_table {"=ACMGRSVTWYHKDBN"};
    return symbol_table[bam_seqi(a_htslib_sequence, index)];
}

HtslibReadFacade::SequenceType HtslibReadFacade::HtslibIterator::get_sequence() const
{
    SequenceType result {};
    auto length = get_sequence_length();
    result.reserve(length);
    auto bam_seq = bam_get_seq(hts_bam1_);
    for (uint32_t i {0}; i < length; ++i) {
        result.push_back(get_base(bam_seq, i));
    }
    return result;
}

std::vector<uint_fast8_t> HtslibReadFacade::HtslibIterator::get_qualities() const
{
    auto qualities = bam_get_qual(hts_bam1_);
    auto length = get_sequence_length();
    std::vector<uint_fast8_t> result {};
    result.reserve(length);
    result.insert(result.begin(), qualities, qualities + length);
    return result;
}

uint32_t HtslibReadFacade::HtslibIterator::get_cigar_length() const noexcept
{
    return hts_bam1_->core.n_cigar;
}

CigarString HtslibReadFacade::HtslibIterator::make_cigar_string() const
{
    auto cigar_operations = bam_get_cigar(hts_bam1_);
    auto length = get_cigar_length();
    
    std::vector<CigarOperation> result;
    result.reserve(length);
    
    for (uint32_t i {0}; i < length; ++i) {
        result.emplace_back(bam_cigar_oplen(cigar_operations[i]), bam_cigar_opchr(cigar_operations[i]));
    };
    
    return CigarString {std::move(result)};
}

std::string HtslibReadFacade::HtslibIterator::get_read_name() const
{
    return std::string {bam_get_qname(hts_bam1_)};
}

std::string HtslibReadFacade::HtslibIterator::get_read_group() const
{
    const auto ptr = bam_aux_get(hts_bam1_.get(), Read_group_tag);
    
    if (ptr == nullptr) {
        throw std::runtime_error {"could not get read group from read " + get_read_name() +
                                    " in file " + hts_facade_.file_path_.string()};
    }
    
    return std::string {bam_aux2Z(ptr)};
}

std::string HtslibReadFacade::HtslibIterator::get_contig_name(HtslibReadFacade::HtsTidType hts_tid) const
{
    return hts_facade_.contig_name_map_.at(hts_tid);
}

// Some of these flags will need to be changes when htslib catches up to the new SAM spec
AlignedRead::FlagData HtslibReadFacade::HtslibIterator::get_flags() const
{
    auto c = hts_bam1_->core;
    AlignedRead::FlagData result {};
    
    result.is_marked_multiple_read_template       = (c.flag & BAM_FPAIRED)        != 0;
    result.is_marked_all_segments_in_read_aligned = (c.flag & BAM_FPROPER_PAIR)   != 0;
    result.is_marked_unmapped                     = (c.flag & BAM_FUNMAP)         != 0;
    result.is_marked_reverse_mapped               = (c.flag & BAM_FREVERSE)       != 0;
    result.is_marked_first_template_segment       = (c.flag & BAM_FREAD1)         != 0;
    result.is_marked_last_template_segmenet       = (c.flag & BAM_FREAD2)         != 0;
    result.is_marked_secondary_alignment          = (c.flag & BAM_FSECONDARY)     != 0;
    result.is_marked_qc_fail                      = (c.flag & BAM_FQCFAIL)        != 0;
    result.is_marked_duplicate                    = (c.flag & BAM_FDUP)           != 0;
    result.is_marked_supplementary_alignment      = (c.flag & BAM_FSUPPLEMENTARY) != 0;
    
    return result;
}

AlignedRead::NextSegment::FlagData HtslibReadFacade::HtslibIterator::get_next_segment_flags() const
{
    auto c = hts_bam1_->core;
    AlignedRead::NextSegment::FlagData result {};
    
    result.is_marked_unmapped       = (c.flag & BAM_FMUNMAP)   != 0;
    result.is_marked_reverse_mapped = (c.flag & BAM_FMREVERSE) != 0;
    
    return result;
}
