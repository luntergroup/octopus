//
//  htslib_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_sam_facade.hpp"

#include <sstream>
#include <cmath>    // std::abs
#include <utility>  // std::move
#include <iterator> // std::begin, std::end, std::next

#include <iostream> // TEST

#include "cigar_string.hpp"

class InvalidBamHeader : std::runtime_error {
public:
    InvalidBamHeader(fs::path file_path, std::string message)
    :
    runtime_error {"invalid BAM header"},
    file_path_ {file_path.string()},
    message_ {std::move(message)}
    {}
    
    const char* what() const noexcept
    {
        return (std::string{runtime_error::what()} + ": in " + file_path_ + " - " + message_).c_str();
    }
    
private:
    std::string message_, file_path_;
};

class InvalidBamRecord : std::runtime_error {
public:
    InvalidBamRecord(fs::path file_path, std::string read_name, std::string message)
    :
    runtime_error {"invalid BAM record"},
    file_path_ {file_path.string()},
    read_name_ {std::move(read_name)},
    message_ {std::move(message)}
    {}
    
    const char* what() const noexcept
    {
        return (std::string {runtime_error::what()} + ": in " + file_path_ + ", read " + read_name_ +
                    " - " + message_).c_str();
    }
    
private:
    std::string message_, read_name_, file_path_;
};

HtslibSamFacade::HtslibSamFacade(const fs::path& file_path)
:
file_path_ {file_path},
hts_file_ {sam_open(file_path_.string().c_str(), "r"), htslib_file_deleter},
hts_header_ {sam_hdr_read(hts_file_.get()), htslib_header_deleter},
hts_index_ {sam_index_load(hts_file_.get(), file_path_.string().c_str()), htslib_index_deleter},
hts_tid_map_ {},
contig_name_map_ {},
sample_map_ {}
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
    
    init_maps();
}

void HtslibSamFacade::open()
{
    hts_file_.reset(sam_open(file_path_.string().c_str(), "r"));
}

void HtslibSamFacade::close()
{
    hts_file_.reset(nullptr);
}

unsigned HtslibSamFacade::get_num_reference_contigs() noexcept
{
    return hts_header_->n_targets;
}

HtslibSamFacade::SizeType HtslibSamFacade::get_reference_contig_size(const std::string& contig_name)
{
    return hts_header_->target_len[get_htslib_tid(contig_name)];
}

uint64_t HtslibSamFacade::get_num_mapped_reads(const std::string& contig_name) const
{
    uint64_t num_mapped, num_unmapped;
    hts_idx_get_stat(hts_index_.get(), get_htslib_tid(contig_name), &num_mapped, &num_unmapped);
    return num_mapped;
}

std::vector<HtslibSamFacade::SampleIdType> HtslibSamFacade::get_samples()
{
    std::vector<HtslibSamFacade::SampleIdType> result {};
    
    for (const auto pair : sample_map_) {
        if (std::find(std::cbegin(result), std::cend(result), pair.second) == std::cend(result)) {
            result.emplace_back(pair.second);
        }
    }
    
    return result;
}

std::vector<std::string> HtslibSamFacade::get_read_groups_in_sample(const SampleIdType& sample)
{
    std::vector<std::string> result {};
    
    for (const auto pair : sample_map_) {
        if (pair.second == sample) result.emplace_back(pair.first);
    }
    
    return result;
}

size_t HtslibSamFacade::count_reads(const GenomicRegion& region)
{
    size_t result {0};
    HtslibIterator it {*this, region};
    
    while (++it) ++result;
    
    return result;
}

size_t HtslibSamFacade::count_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    size_t result {0};
    HtslibIterator it {*this, region};
    
    while (++it && sample_map_.at(it.get_read_group()) == sample) ++result;
    
    return result;
}

GenomicRegion HtslibSamFacade::find_head_region(const GenomicRegion& region, size_t target_coverage)
{
    HtslibIterator it {*this, region};
    
    while (++it && target_coverage > 0) --target_coverage;
    
    return GenomicRegion {region.get_contig_name(), region.get_begin(), (*it).first.get_region().get_end()};
}

HtslibSamFacade::SampleIdToReadsMap HtslibSamFacade::fetch_reads(const GenomicRegion& region)
{
    HtslibIterator it {*this, region};
    SampleIdToReadsMap result {};
    
    while (++it) {
        try {
            auto a_read_and_its_group = *it;
            const auto& sample = sample_map_.at(a_read_and_its_group.second);
            result[std::move(sample)].emplace_back(std::move(a_read_and_its_group.first));
        } catch (InvalidBamRecord& e) {
            // TODO: Just ignore? Could log or something.
            //std::clog << "Warning: " << e.what() << std::endl;
        } catch (...) {
            // Could be something really bad.
            throw;
        }
    }
    
    for (auto& sample_reads_pair : result) {
        sample_reads_pair.second.shrink_to_fit();
    }
    
    return result;
}

std::vector<AlignedRead> HtslibSamFacade::fetch_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    HtslibIterator it {*this, region};
    std::vector<AlignedRead> result {};
    
    while (++it) {
        try {
            if (sample_map_.at(it.get_read_group()) == sample) {
                result.push_back((*it).first);
            }
        } catch (InvalidBamRecord& e) {
            // TODO: Just ignore? Could log or something.
            //std::clog << "Warning: " << e.what() << std::endl;
        } catch (...) {
            // Could be something really bad.
            throw;
        }
    }
    
    return result;
}

std::vector<std::string> HtslibSamFacade::get_reference_contig_names()
{
    std::vector<std::string> result {};
    result.reserve(get_num_reference_contigs());
    
    for (HtsTidType hts_tid {0}; hts_tid < get_num_reference_contigs(); ++hts_tid) {
        result.emplace_back(get_contig_name(hts_tid));
    }
    
    return result;
}

std::vector<GenomicRegion> HtslibSamFacade::get_possible_regions_in_file()
{
    std::vector<GenomicRegion> result {};
    result.reserve(get_num_reference_contigs());
    
    for (HtsTidType hts_tid {0}; hts_tid < get_num_reference_contigs(); ++hts_tid) {
        auto contig_name = get_contig_name(hts_tid);
        // CRAM files don't seem to have the same index stats as BAM files so
        // we don't know which contigs have been mapped to
        if (hts_file_->is_cram || get_num_mapped_reads(contig_name) > 0) {
            result.emplace_back(std::move(contig_name), 0, get_reference_contig_size(contig_name));
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

bool is_tag_type(const std::string& header_line, const char* tag)
{
    return header_line.compare(1, 2, tag) == 0;
}

bool has_tag(const std::string& header_line, const char* tag)
{
    return header_line.find(tag) != std::string::npos;
}

std::string get_tag_value(const std::string& line, const char* tag)
{
    // format TAG:VALUE\t
    auto tag_position = line.find(tag);
    
    if (tag_position != std::string::npos) {
        auto value_position = line.find(':', tag_position) + 1;
        auto tag_value_size = line.find('\t', value_position) - value_position;
        return line.substr(value_position, tag_value_size);
    } else {
        throw std::runtime_error {"no " + std::string {tag} + " tag"};
    }
}

void HtslibSamFacade::init_maps()
{
    hts_tid_map_.reserve(get_num_reference_contigs());
    contig_name_map_.reserve(get_num_reference_contigs());
    
    for (HtsTidType hts_tid {0}; hts_tid < get_num_reference_contigs(); ++hts_tid) {
        hts_tid_map_.emplace(hts_header_->target_name[hts_tid], hts_tid);
        contig_name_map_.emplace(hts_tid, hts_header_->target_name[hts_tid]);
    }
    
    std::string the_header_text (hts_header_->text, hts_header_->l_text);
    std::stringstream ss {the_header_text};
    std::string line {};
    unsigned num_read_groups {0};
    
    while (std::getline(ss, line, '\n')) {
        if (is_tag_type(line, Read_group_tag)) {
            if (!has_tag(line, Read_group_id_tag)) {
                throw InvalidBamHeader {file_path_, "no read group identifier tag (ID) in @RG line"};
            }
            if (!has_tag(line, Sample_id_tag)) {
                // The SAM specification does not specify the sample tag 'SM' as a required,
                // however we can't do much without it.
                throw InvalidBamHeader {file_path_, "no sample tag (SM) in @RG line"};
            }
            sample_map_.emplace(get_tag_value(line, Read_group_id_tag), get_tag_value(line, Sample_id_tag));
            ++num_read_groups;
        }
    }
    
    if (num_read_groups == 0) {
        throw InvalidBamHeader {file_path_, "no read group (@RG) lines found"};
    }
}

HtslibSamFacade::HtsTidType HtslibSamFacade::get_htslib_tid(const std::string& contig_name) const
{
    return hts_tid_map_.at(contig_name);
}

const std::string& HtslibSamFacade::get_contig_name(HtsTidType hts_tid) const
{
    return contig_name_map_.at(hts_tid);
}

// HtslibIterator

HtslibSamFacade::HtslibIterator::HtslibIterator(HtslibSamFacade& hts_facade, const GenomicRegion& region)
:
hts_facade_ {hts_facade},
hts_iterator_ {sam_itr_querys(hts_facade_.hts_index_.get(), hts_facade_.hts_header_.get(),
                              to_string(region).c_str()), htslib_iterator_deleter},
hts_bam1_ {bam_init1(), htslib_bam1_deleter}
{
    if (hts_iterator_ == nullptr) {
        throw std::runtime_error {"could not load read iterator for " + hts_facade.file_path_.string()};
    }
    
    if (hts_bam1_ == nullptr) {
        throw std::runtime_error {"error creating bam1 for " + hts_facade.file_path_.string()};
    }
}

std::string get_read_name(bam1_t* b)
{
    return std::string {bam_get_qname(b)};
}

HtslibSamFacade::ReadGroupIdType HtslibSamFacade::HtslibIterator::get_read_group() const
{
    const auto ptr = bam_aux_get(hts_bam1_.get(), Read_group_tag);
    
    if (ptr == nullptr) {
        throw InvalidBamRecord {hts_facade_.file_path_, get_read_name(hts_bam1_.get()), "no read group"};
    }
    
    return HtslibSamFacade::ReadGroupIdType {bam_aux2Z(ptr)};
}

bool HtslibSamFacade::HtslibIterator::operator++()
{
    return sam_itr_next(hts_facade_.hts_file_.get(), hts_iterator_.get(), hts_bam1_.get()) >= 0;
}

auto get_read_pos(bam1_t* b) noexcept
{
    return b->core.pos;
}

auto get_sequence_length(bam1_t* b) noexcept
{
    return b->core.l_qseq;
}

char get_base(uint8_t* hts_sequence, uint32_t index) noexcept
{
    static constexpr const char* symbol_table {"=ACMGRSVTWYHKDBN"};
    return symbol_table[bam_seqi(hts_sequence, index)];
}

AlignedRead::SequenceType get_sequence(bam1_t* b)
{
    AlignedRead::SequenceType result {};
    auto length = get_sequence_length(b);
    result.reserve(length);
    auto bam_seq = bam_get_seq(b);
    for (uint32_t i {0}; i < length; ++i) {
        result.push_back(get_base(bam_seq, i));
    }
    return result;
}

std::vector<AlignedRead::QualityType> get_qualities(bam1_t* b)
{
    auto qualities = bam_get_qual(b);
    auto length = get_sequence_length(b);
    std::vector<AlignedRead::QualityType> result {};
    result.reserve(length);
    result.insert(result.begin(), qualities, qualities + length);
    return result;
}

auto get_cigar_length(bam1_t* b) noexcept
{
    return b->core.n_cigar;
}

CigarString get_cigar_string(bam1_t* b)
{
    auto cigar_operations = bam_get_cigar(b);
    auto length = get_cigar_length(b);
    
    std::vector<CigarOperation> result;
    result.reserve(length);
    
    for (uint32_t i {0}; i < length; ++i) {
        result.emplace_back(bam_cigar_oplen(cigar_operations[i]), bam_cigar_opchr(cigar_operations[i]));
    };
    
    return CigarString {std::move(result)};
}

// Some of these flags will need to be changes when htslib catches up to the new SAM spec
AlignedRead::FlagData get_flags(bam1_t* b)
{
    auto c = b->core;
    
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

AlignedRead::NextSegment::FlagData get_next_segment_flags(bam1_t* b)
{
    auto c = b->core;
    AlignedRead::NextSegment::FlagData result {};
    
    result.is_marked_unmapped       = (c.flag & BAM_FMUNMAP)   != 0;
    result.is_marked_reverse_mapped = (c.flag & BAM_FMREVERSE) != 0;
    
    return result;
}

std::pair<AlignedRead, HtslibSamFacade::ReadGroupIdType> HtslibSamFacade::HtslibIterator::operator*() const
{
    auto qualities = get_qualities(hts_bam1_.get());
    
    if (qualities.empty() || qualities[0] == 0xff) {
        throw InvalidBamRecord {hts_facade_.file_path_, get_read_name(hts_bam1_.get()), "corrupt sequence data"};
    }
    
    auto cigar = get_cigar_string(hts_bam1_.get());
    
    if (cigar.empty()) {
        throw InvalidBamRecord {hts_facade_.file_path_, get_read_name(hts_bam1_.get()), "empty cigar string"};
    }
    
    auto c = hts_bam1_->core;
    
    auto read_begin_tmp = soft_clipped_read_begin(cigar, c.pos);
    
    auto sequence = get_sequence(hts_bam1_.get());
    
    if (read_begin_tmp < 0) {
        // i.e. if the read hangs off the left of the contig
        auto overhang = std::abs(read_begin_tmp);
        sequence.erase(std::begin(sequence), std::next(std::begin(sequence), overhang));
        qualities.erase(std::begin(qualities), std::next(std::begin(qualities), overhang));
        read_begin_tmp = 0;
    }
    
    auto read_begin = static_cast<AlignedRead::SizeType>(read_begin_tmp);
    
    const auto& contig_name = hts_facade_.get_contig_name(c.tid);
    
    if (c.mtid == -1) { // i.e. has no mate TODO: check if this is always true
        return {AlignedRead {
            GenomicRegion {contig_name, read_begin, read_begin + reference_size<GenomicRegion::SizeType>(cigar)},
            std::move(sequence),
            std::move(qualities),
            std::move(cigar),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags(hts_bam1_.get())
        }, get_read_group()};
    } else {
        return {AlignedRead {
            GenomicRegion {contig_name, read_begin, read_begin + reference_size<GenomicRegion::SizeType>(cigar)},
            std::move(sequence),
            std::move(qualities),
            std::move(cigar),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags(hts_bam1_.get()),
            contig_name,
            static_cast<AlignedRead::SizeType>(c.mpos),
            static_cast<AlignedRead::SizeType>(std::abs(c.isize)),
            get_next_segment_flags(hts_bam1_.get())
        }, get_read_group()};
    }
}
