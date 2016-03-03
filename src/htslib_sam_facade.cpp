//
//  htslib_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_sam_facade.hpp"

#include <sstream>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <stdexcept>

#include "cigar_string.hpp"
#include "genomic_region.hpp"
#include "contig_region.hpp"

#include <iostream> // TEST

class InvalidBamHeader : std::runtime_error {
public:
    InvalidBamHeader(boost::filesystem::path file_path, std::string message)
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
    InvalidBamRecord(boost::filesystem::path file_path, std::string read_name, std::string message)
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

// public methods

HtslibSamFacade::HtslibSamFacade(Path file_path)
:
file_path_ {std::move(file_path)},
hts_file_ {sam_open(file_path_.c_str(), "r"), HtsFileDeleter {}},
hts_header_ {(hts_file_) ? sam_hdr_read(hts_file_.get()) : nullptr, HtsHeaderDeleter {}},
hts_index_ {(hts_file_) ? sam_index_load(hts_file_.get(), file_path_.c_str()) : nullptr, HtsIndexDeleter {}},
hts_tids_ {},
contig_names_ {},
sample_names_ {},
samples_ {}
{
    if (is_open()) {
        init_maps();
        
        for (const auto& pair : sample_names_) {
            if (std::find(std::cbegin(samples_), std::cend(samples_), pair.second) == std::cend(samples_)) {
                samples_.emplace_back(pair.second);
            }
        }
        
        samples_.shrink_to_fit();
    } else {
        close();
    }
}

bool HtslibSamFacade::is_open() const noexcept
{
    return hts_file_ != nullptr && hts_header_ != nullptr && hts_index_ != nullptr;
}

void HtslibSamFacade::open()
{
    hts_file_.reset(sam_open(file_path_.string().c_str(), "r"));
    if (hts_file_) {
        hts_header_.reset(sam_hdr_read(hts_file_.get()));
        hts_index_.reset(sam_index_load(hts_file_.get(), file_path_.c_str()));
    }
}

void HtslibSamFacade::close()
{
    hts_file_.reset(nullptr);
    hts_header_.reset(nullptr);
    hts_index_.reset(nullptr);
}

unsigned HtslibSamFacade::count_reference_contigs()
{
    if (!is_open()) {
        throw std::runtime_error {"HtslibSamFacade: file not open"};
    }
    return hts_header_->n_targets;
}

HtslibSamFacade::SizeType HtslibSamFacade::get_reference_contig_size(const std::string& contig_name)
{
    return hts_header_->target_len[get_htslib_tid(contig_name)];
}

uint64_t HtslibSamFacade::get_num_mapped_reads(const std::string& contig_name) const
{
    uint64_t num_mapped {}, num_unmapped {};
    hts_idx_get_stat(hts_index_.get(), get_htslib_tid(contig_name), &num_mapped, &num_unmapped);
    return num_mapped;
}

std::vector<HtslibSamFacade::SampleIdType> HtslibSamFacade::extract_samples()
{
    return samples_;
}

std::vector<std::string> HtslibSamFacade::extract_read_groups_in_sample(const SampleIdType& sample)
{
    std::vector<std::string> result {};
    
    for (const auto pair : sample_names_) {
        if (pair.second == sample) result.emplace_back(pair.first);
    }
    
    result.shrink_to_fit();
    
    return result;
}

bool HtslibSamFacade::has_contig_reads(const GenomicRegion::ContigNameType& contig)
{
    HtslibIterator it {*this, contig};
    return ++it;
}

size_t HtslibSamFacade::count_reads(const GenomicRegion& region)
{
    HtslibIterator it {*this, region};
    
    size_t result {0};
    
    while (++it) ++result;
    
    return result;
}

size_t HtslibSamFacade::count_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    HtslibIterator it {*this, region};
    
    size_t result {0};
    
    while (++it && sample_names_.at(it.get_read_group()) == sample) ++result;
    
    return result;
}

std::pair<GenomicRegion, std::vector<unsigned>>
HtslibSamFacade::find_covered_subregion(const GenomicRegion& region, size_t max_coverage)
{
    HtslibIterator it {*this, region};
    
    if (max_coverage == 0 || !++it) {
        return std::make_pair(head_region(region), std::vector<unsigned> {});
    }
    
    --max_coverage;
    
    const auto first_begin = it.get_begin();
    
    std::vector<unsigned> position_counts(max_coverage, 0);
    
    ++position_counts[0];
    
    auto read_begin = first_begin;
    
    while (++it && --max_coverage > 0) {
        read_begin = it.get_begin();
        
        const auto offset = read_begin - first_begin;
        
        if (offset >= position_counts.size()) {
            position_counts.resize(std::max(offset, position_counts.size() + max_coverage));
        } else {
            ++position_counts[offset];
        }
    }
    
    position_counts.resize(read_begin - first_begin);
    
    const auto result_begin = static_cast<GenomicRegion::SizeType>(first_begin);
    const auto result_end   = static_cast<GenomicRegion::SizeType>(read_begin);
    
    return std::make_pair(GenomicRegion {region.get_contig_name(), result_begin, result_end},
                          position_counts);
}

bool contains(const std::vector<HtslibSamFacade::SampleIdType>& samples,
              const HtslibSamFacade::SampleIdType& sample)
{
    return std::find(std::cbegin(samples), std::cend(samples), sample) != std::cend(samples);
}

bool is_subset(std::vector<HtslibSamFacade::SampleIdType> lhs,
               std::vector<HtslibSamFacade::SampleIdType> rhs)
{
    std::sort(std::begin(lhs), std::end(lhs));
    std::sort(std::begin(rhs), std::end(rhs));
    return std::includes(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::cend(rhs));
}

std::pair<GenomicRegion, std::vector<unsigned>>
HtslibSamFacade::find_covered_subregion(const std::vector<SampleIdType>& samples,
                                        const GenomicRegion& region, size_t max_coverage)
{
    if (is_subset(samples, samples_)) {
        return find_covered_subregion(region, max_coverage);
    }
    
    HtslibIterator it {*this, region};
    
    if (max_coverage == 0 || !++it) {
        return std::make_pair(head_region(region), std::vector<unsigned> {});
    }
    
    while (!contains(samples, sample_names_.at(it.get_read_group()))) {
        if (!++it) {
            return std::make_pair(head_region(region), std::vector<unsigned> {});
        }
    };
    
    --max_coverage;
    
    const auto first_begin = it.get_begin();
    
    std::vector<unsigned> position_counts(max_coverage, 0);
    
    ++position_counts[0];
    
    auto read_begin = first_begin;
    
    while (++it && !contains(samples, sample_names_.at(it.get_read_group())) && --max_coverage > 0) {
        read_begin = it.get_begin();
        
        const auto offset = read_begin - first_begin;
        
        if (offset >= position_counts.size()) {
            position_counts.resize(std::max(offset, position_counts.size() + max_coverage));
        } else {
            ++position_counts[offset];
        }
    }
    
    position_counts.resize(read_begin - first_begin);
    
    const auto result_begin = static_cast<GenomicRegion::SizeType>(first_begin);
    const auto result_end   = static_cast<GenomicRegion::SizeType>(read_begin);
    
    return std::make_pair(GenomicRegion {region.get_contig_name(), result_begin, result_end},
                          position_counts);
}

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const GenomicRegion& region)
{
    SampleReadMap result {samples_.size()};
    
    if (samples_.size() == 1) {
        auto reads = fetch_reads(samples_.front(), region);
        result.emplace(samples_.front(), std::move(reads));
        return result;
    }
    
    HtslibIterator it {*this, region};
    
    for (const auto& sample : samples_) {
        auto p = result.emplace(std::piecewise_construct, std::forward_as_tuple(sample),
                                std::forward_as_tuple());
        p.first->second.reserve(default_reserve_);
    }
    
    while (++it) {
        auto read = *it;
        if (read) {
             result[sample_names_.at(it.get_read_group())].emplace(*std::move(read));
        }
    }
    
    for (auto& p : result) p.second.shrink_to_fit();
    
    return result;
}

HtslibSamFacade::Reads HtslibSamFacade::fetch_reads(const SampleIdType& sample, const GenomicRegion& region)
{
    HtslibIterator it {*this, region};
    
    Reads result {};
    
    if (std::find(std::cbegin(samples_), std::cend(samples_), sample) == std::cend(samples_)) {
        return result;
    }
    
    result.reserve(default_reserve_);
    
    if (samples_.size() == 1) {
        while (++it) {
            auto read = *it;
            if (read) {
                result.emplace(*std::move(read));
            }
        }
    } else {
        while (++it) {
            if (sample_names_.at(it.get_read_group()) == sample) {
                auto read = *it;
                if (read) {
                    result.emplace(*std::move(read));
                }
            }
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

HtslibSamFacade::SampleReadMap HtslibSamFacade::fetch_reads(const std::vector<SampleIdType>& samples,
                                                            const GenomicRegion& region)
{
    if (is_subset(samples_, samples)) return fetch_reads(region);
    
    HtslibIterator it {*this, region};
    
    SampleReadMap result {samples.size()};
    
    for (const auto& sample : samples) {
        if (std::find(std::cbegin(samples_), std::cend(samples_), sample) != std::cend(samples_)) {
            auto p = result.emplace(std::piecewise_construct, std::forward_as_tuple(sample),
                                    std::forward_as_tuple());
            p.first->second.reserve(default_reserve_);
        }
    }
    
    if (result.empty()) return result; // no matching samples
    
    while (++it) {
        const auto& sample = sample_names_.at(it.get_read_group());
        if (result.count(sample) == 1) {
            auto read = *it;
            if (read) {
                result[sample].emplace(*std::move(read));
            }
        }
    }
    
    for (auto& p : result) p.second.shrink_to_fit();
    
    return result;
}

std::vector<std::string> HtslibSamFacade::extract_reference_contig_names()
{
    std::vector<std::string> result {};
    result.reserve(count_reference_contigs());
    
    for (HtsTidType hts_tid {}; hts_tid < count_reference_contigs(); ++hts_tid) {
        result.emplace_back(get_contig_name(hts_tid));
    }
    
    return result;
}

std::vector<GenomicRegion> HtslibSamFacade::extract_possible_regions_in_file()
{
    std::vector<GenomicRegion> result {};
    result.reserve(count_reference_contigs());
    
    for (HtsTidType hts_tid {0}; hts_tid < count_reference_contigs(); ++hts_tid) {
        const auto& contig_name = get_contig_name(hts_tid);
        // CRAM files don't seem to have the same index stats as BAM files so
        // we don't know which contigs have been mapped to
        if (hts_file_->is_cram || get_num_mapped_reads(contig_name) > 0) {
            // without actually looking at the reads we don't know what coverage there is so
            // we just assume the entire contig is possible
            result.emplace_back(contig_name, 0, get_reference_contig_size(contig_name));
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

// private methods

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
    hts_tids_.reserve(count_reference_contigs());
    contig_names_.reserve(count_reference_contigs());
    
    for (HtsTidType hts_tid {}; hts_tid < count_reference_contigs(); ++hts_tid) {
        hts_tids_.emplace(hts_header_->target_name[hts_tid], hts_tid);
        contig_names_.emplace(hts_tid, hts_header_->target_name[hts_tid]);
    }
    
    std::string the_header_text (hts_header_->text, hts_header_->l_text);
    std::stringstream ss {the_header_text};
    std::string line {};
    unsigned num_read_groups {};
    
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
            
            sample_names_.emplace(get_tag_value(line, Read_group_id_tag), get_tag_value(line, Sample_id_tag));
            ++num_read_groups;
        }
    }
    
    if (num_read_groups == 0) {
        throw InvalidBamHeader {file_path_, "no read group (@RG) lines found"};
    }
}

HtslibSamFacade::HtsTidType HtslibSamFacade::get_htslib_tid(const std::string& contig_name) const
{
    return hts_tids_.at(contig_name);
}

const std::string& HtslibSamFacade::get_contig_name(HtsTidType hts_tid) const
{
    return contig_names_.at(hts_tid);
}

// HtslibIterator

HtslibSamFacade::HtslibIterator::HtslibIterator(HtslibSamFacade& hts_facade, const GenomicRegion& region)
:
hts_facade_ {hts_facade},
hts_iterator_ {hts_facade.is_open() ? sam_itr_querys(hts_facade_.hts_index_.get(), hts_facade_.hts_header_.get(),
                                                     to_string(region).c_str()) : nullptr, HtsIteratorDeleter {}},
hts_bam1_ {bam_init1(), HtsBam1Deleter {}}
{
    if (hts_iterator_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: could not load iterator for " + hts_facade.file_path_.string()};
    }
    
    if (hts_bam1_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: error creating bam1 for " + hts_facade.file_path_.string()};
    }
}

HtslibSamFacade::HtslibIterator::HtslibIterator(HtslibSamFacade& hts_facade, const GenomicRegion::ContigNameType& contig)
:
hts_facade_ {hts_facade},
hts_iterator_ {hts_facade.is_open() ? sam_itr_querys(hts_facade_.hts_index_.get(), hts_facade_.hts_header_.get(),
                                                     contig.c_str()) : nullptr, HtsIteratorDeleter {}},
hts_bam1_ {bam_init1(), HtsBam1Deleter {}}
{
    if (hts_iterator_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: could not load iterator for " + hts_facade.file_path_.string()};
    }
    
    if (hts_bam1_ == nullptr) {
        throw std::runtime_error {"HtslibIterator: error creating bam1 for " + hts_facade.file_path_.string()};
    }
}

std::string get_read_name(const bam1_t* b)
{
    return std::string {bam_get_qname(b)};
}

bool HtslibSamFacade::HtslibIterator::operator++()
{
    return sam_itr_next(hts_facade_.hts_file_.get(), hts_iterator_.get(), hts_bam1_.get()) >= 0;
}

auto get_read_pos(const bam1_t* b) noexcept
{
    return b->core.pos;
}

auto get_sequence_length(const bam1_t* b) noexcept
{
    return b->core.l_qseq;
}

char get_base(const uint8_t* hts_sequence, const uint32_t index) noexcept
{
    static constexpr const char* symbol_table {"=ACMGRSVTWYHKDBN"};
    return symbol_table[bam_seqi(hts_sequence, index)];
}

AlignedRead::SequenceType get_sequence(const bam1_t* b)
{
    using SequenceType = AlignedRead::SequenceType;
    
    const auto sequence_length  = static_cast<SequenceType::size_type>(get_sequence_length(b));
    const auto hts_sequence     = bam_get_seq(b);
    
    SequenceType result(sequence_length, 'N');
    
    uint32_t i {};
    std::generate_n(std::begin(result), sequence_length,
                    [&i, &hts_sequence] () { return get_base(hts_sequence, i++); });
    
    return result;
}

AlignedRead::Qualities get_qualities(const bam1_t* b)
{
    const auto qualities = bam_get_qual(b);
    const auto length    = get_sequence_length(b);
    return AlignedRead::Qualities(qualities, qualities + length);
}

auto get_cigar_length(const bam1_t* b) noexcept
{
    return b->core.n_cigar;
}

CigarString get_cigar_string(const bam1_t* b)
{
    const auto cigar_operations = bam_get_cigar(b);
    const auto cigar_length     = get_cigar_length(b);
    
    CigarString result(cigar_length);
    
    std::transform(cigar_operations, cigar_operations + cigar_length, std::begin(result),
                   [] (auto op) {
                       return CigarOperation {bam_cigar_oplen(op), bam_cigar_opchr(op)};
                   });
    
    return result;
}

// Some of these flags will need to be changes when htslib catches up to the new SAM spec
AlignedRead::Flags get_flags(const bam1_t* b) noexcept
{
    const auto c = b->core;
    
    AlignedRead::Flags result {};
    
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

AlignedRead::NextSegment::Flags get_next_segment_flags(const bam1_t* b)
{
    const auto c = b->core;
    
    AlignedRead::NextSegment::Flags result {};
    
    result.is_marked_unmapped       = (c.flag & BAM_FMUNMAP)   != 0;
    result.is_marked_reverse_mapped = (c.flag & BAM_FMREVERSE) != 0;
    
    return result;
}

boost::optional<AlignedRead> HtslibSamFacade::HtslibIterator::operator*() const
{
    using std::begin; using std::end; using std::next; using std::move;
    
    auto qualities = get_qualities(hts_bam1_.get());
    
    if (qualities.empty() || qualities[0] == 0xff) return boost::none;
    
    auto cigar = get_cigar_string(hts_bam1_.get());
    
    if (!is_valid_cigar(cigar)) return boost::none;
    
    const auto c = hts_bam1_->core;
    
    auto read_begin_tmp = soft_clipped_read_begin(cigar, c.pos);
    
    auto sequence = get_sequence(hts_bam1_.get());
    
    if (read_begin_tmp < 0) {
        // Then the read hangs off the left of the contig, and we must remove bases, qualities, and
        // adjust the cigar string as we cannot have a negative begin position
        
        auto overhang_size = std::abs(read_begin_tmp);
        
        sequence.erase(begin(sequence), next(begin(sequence), overhang_size));
        qualities.erase(begin(qualities), next(begin(qualities), overhang_size));
        
        auto soft_clip_size = cigar.front().get_size();
        
        if (overhang_size == soft_clip_size) {
            cigar.erase(begin(cigar));
        } else { // then soft_clip_size > overhang_size
            cigar.front() = CigarOperation {soft_clip_size - overhang_size, CigarOperation::SOFT_CLIPPED};
        }
        
        read_begin_tmp = 0;
    }
    
    const auto read_begin = static_cast<AlignedRead::SizeType>(read_begin_tmp);
    
    const auto& contig_name = hts_facade_.get_contig_name(c.tid);
    
    if (c.mtid == -1) { // i.e. has no mate TODO: check if this is always true
        return AlignedRead {
            GenomicRegion {contig_name, read_begin, read_begin + reference_size<GenomicRegion::SizeType>(cigar)},
            move(sequence),
            move(qualities),
            move(cigar),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags(hts_bam1_.get())
        };
    } else {
        return AlignedRead {
            GenomicRegion {contig_name, read_begin, read_begin + reference_size<GenomicRegion::SizeType>(cigar)},
            move(sequence),
            move(qualities),
            move(cigar),
            static_cast<AlignedRead::QualityType>(c.qual),
            get_flags(hts_bam1_.get()),
            contig_name,
            static_cast<AlignedRead::SizeType>(c.mpos),
            static_cast<AlignedRead::SizeType>(std::abs(c.isize)),
            get_next_segment_flags(hts_bam1_.get())
        };
    }
}

HtslibSamFacade::ReadGroupIdType HtslibSamFacade::HtslibIterator::get_read_group() const
{
    const auto ptr = bam_aux_get(hts_bam1_.get(), Read_group_tag);
    
    if (ptr == nullptr) {
        throw InvalidBamRecord {hts_facade_.file_path_, get_read_name(hts_bam1_.get()), "no read group"};
    }
    
    return HtslibSamFacade::ReadGroupIdType {bam_aux2Z(ptr)};
}

bool HtslibSamFacade::HtslibIterator::is_good() const noexcept
{
    if (get_sequence_length(hts_bam1_.get()) == 0) {
        return false;
    }
    
    if (bam_get_qual(hts_bam1_.get())[0] == 0xff) {
        return false;
    }
    
    const auto cigar_length = get_cigar_length(hts_bam1_.get());
    
    if (cigar_length == 0) {
        return false;
    }
    
    const auto cigar_operations = bam_get_cigar(hts_bam1_.get());
    
    return std::all_of(cigar_operations, cigar_operations + cigar_length,
                       [] (const auto op) {
                           return bam_cigar_oplen(op) > 0;
                       });
}

std::size_t HtslibSamFacade::HtslibIterator::get_begin() const noexcept
{
    return static_cast<ContigRegion::SizeType>(hts_bam1_->core.pos);
}
