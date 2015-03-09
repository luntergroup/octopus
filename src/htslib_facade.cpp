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

uint_fast32_t HtslibFacade::get_num_reference_contigs() noexcept
{
    return the_header_->n_targets;
}

uint_fast32_t HtslibFacade::get_reference_contig_size(const std::string& contig_name)
{
    return the_header_->target_len[get_htslib_tid(contig_name)];
}

int HtslibFacade::HtslibIterator::operator++()
{
    return sam_itr_next(hts_facade_.the_file_.get(), the_iterator_.get(), the_bam1_.get()) >= 0;
}

uint_fast32_t HtslibFacade::HtslibIterator::get_read_start() const noexcept
{
    return static_cast<uint_fast32_t>(the_bam1_->core.pos);
}

uint32_t HtslibFacade::HtslibIterator::get_sequence_length() const noexcept
{
    return the_bam1_->core.l_qseq;
}

char HtslibFacade::HtslibIterator::get_base(uint8_t* a_htslib_sequence,
                                                   uint32_t index) const noexcept
{
    static constexpr const char* symbol_table {"=ACMGRSVTWYHKDBN"};
    return symbol_table[bam_seqi(a_htslib_sequence, index)];
}

std::string HtslibFacade::HtslibIterator::get_sequence() const
{
    std::string result {};
    auto length = get_sequence_length();
    result.reserve(length);
    for (uint32_t i = 0; i < length; ++i) {
        result.push_back(get_base(bam_get_seq(the_bam1_), i));
    }
    return result;
}

std::vector<uint_fast8_t> HtslibFacade::HtslibIterator::get_qualities() const
{
    auto qualities = bam_get_qual(the_bam1_);
    auto length = get_sequence_length();
    std::vector<uint_fast8_t> result;
    result.reserve(length);
    result.insert(result.begin(), qualities, qualities + length);
    return result;
}

uint32_t HtslibFacade::HtslibIterator::get_cigar_length() const noexcept
{
    return the_bam1_->core.n_cigar;
}

CigarString HtslibFacade::HtslibIterator::make_cigar_string() const
{
    static constexpr const char* operation_table {"MIDNSHP=X"};
    auto cigar_operations = bam_get_cigar(the_bam1_);
    auto length = get_cigar_length();
    std::vector<CigarOperation> result;
    result.reserve(length);
    for (uint32_t i {0}; i < length; ++i) {
        result.emplace_back(bam_cigar_oplen(cigar_operations[i]),
                            operation_table[bam_cigar_op(cigar_operations[i])]);
    };
    return CigarString(std::move(result));
}

std::string HtslibFacade::HtslibIterator::get_read_group() const
{
    return std::string(bam_aux2Z(bam_aux_get(the_bam1_.get(), Read_group_tag)));
}

std::string HtslibFacade::HtslibIterator::get_contig_name(int32_t htslib_tid) const
{
    return hts_facade_.contig_name_map_.at(htslib_tid);
}

std::string HtslibFacade::HtslibIterator::get_read_name() const
{
    return std::string {bam_get_qname(the_bam1_)};
}

std::string HtslibFacade::get_reference_contig_name(int32_t hts_tid) const
{
    return std::string(the_header_->target_name[hts_tid]);
}

bool HtslibFacade::is_type(const std::string& header_line, const char* tag) const
{
    return header_line.compare(1, 2, tag) == 0;
}

bool HtslibFacade::has_tag(const std::string& header_line, const char* tag) const
{
    return header_line.find(tag) != std::string::npos;
}

std::string HtslibFacade::get_tag_value(const std::string& line, const char* tag) const
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

uint64_t HtslibFacade::get_num_mapped_reads(const std::string& reference_contig_name) const
{
    uint64_t num_mapped, num_unmapped;
    hts_idx_get_stat(the_index_.get(), get_htslib_tid(reference_contig_name), &num_mapped, &num_unmapped);
    return num_mapped;
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

std::size_t HtslibFacade::get_num_reads(const GenomicRegion& a_region)
{
    std::size_t result {0};
    HtslibIterator it {*this, a_region};
    while (++it) { ++result; }
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
    auto the_cigar_string = make_cigar_string();
    if (the_cigar_string.empty()) throw std::runtime_error {"bad sequence"};
    auto c = the_bam1_->core;
    auto read_start = static_cast<AlignedRead::SizeType>(get_soft_clipped_read_begin(the_cigar_string, c.pos));
    
    if (c.mtid == -1) { // TODO: check if this is always true
        return {AlignedRead {
            GenomicRegion(get_contig_name(c.tid), read_start, read_start + get_sequence_length()),
            get_sequence(),
            std::move(the_qualities),
            std::move(the_cigar_string),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags()
        }, get_read_group()};
    } else {
        return {AlignedRead {
            GenomicRegion(get_contig_name(c.tid), read_start, read_start + get_sequence_length()),
            get_sequence(),
            std::move(the_qualities),
            std::move(the_cigar_string),
            static_cast<AlignedRead::MatePair::InsertSizeType>(c.isize),
            get_contig_name(c.mtid),
            static_cast<AlignedRead::SizeType>(c.mpos),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags(),
            get_mate_flags()
        }, get_read_group()};
    }
}

AlignedRead::SupplementaryData HtslibFacade::HtslibIterator::get_flags() const
{
    auto c = the_bam1_->core;
    AlignedRead::SupplementaryData result {};
    result.is_marked_duplicate           = (c.flag & BAM_FDUP)         != 0;
    result.is_marked_unmapped            = (c.flag & BAM_FUNMAP)       != 0;
    result.is_marked_reverse_mapped      = (c.flag & BAM_FREVERSE)     != 0;
    result.is_marked_paired              = (c.flag & BAM_FPAIRED)      != 0;
    result.is_marked_proper_pair         = (c.flag & BAM_FPROPER_PAIR) != 0;
    result.is_marked_secondary_alignment = (c.flag & BAM_FSECONDARY)   != 0;
    result.is_marked_qc_fail             = (c.flag & BAM_FQCFAIL)      != 0;
    return result;
}

AlignedRead::MatePair::SupplementaryData HtslibFacade::HtslibIterator::get_mate_flags() const
{
    auto c = the_bam1_->core;
    AlignedRead::MatePair::SupplementaryData result {};
    result.is_marked_unmapped       = (c.flag & BAM_FMUNMAP)   != 0;
    result.is_marked_reverse_mapped = (c.flag & BAM_FMREVERSE) != 0;
    return result;
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
