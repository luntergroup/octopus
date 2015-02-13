//
//  htslib_facade.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "htslib_facade.h"

#include <stdexcept>

std::set<AlignedRead> HtslibFacade::fetch_reads(const GenomicRegion& a_region)
{
    HtslibIterator it {*this, a_region};
    std::set<AlignedRead> the_reads {};
    
    while (++it) {
        the_reads.emplace(*it);
    }
    
    return the_reads;
}

std::unordered_set<std::string> HtslibFacade::get_reference_contig_names()
{
    std::unordered_set<std::string> result {};
    result.reserve(get_num_reference_contigs());
    for (uint_fast32_t i {0}; i < get_num_reference_contigs(); ++i) {
        result.emplace(the_header_->target_name[i]);
    }
    return result;
}

std::unordered_map<int32_t, std::string>
HtslibFacade::get_htslib_tid_to_contig_name_mappings() const
{
    std::unordered_map<int32_t, std::string> result {};
    result.reserve(the_header_->n_targets);
    for (uint_fast32_t i {0}; i < the_header_->n_targets; ++i) {
        result.emplace(std::make_pair(i, the_header_->target_name[i]));
    }
    return result;
}

AlignedRead HtslibFacade::HtslibIterator::operator*() const
{
    bam1_core_t c = b_->core;
    uint_fast32_t sequence_length = c.l_qseq;
    
    if (sequence_length == 0) {
        // TODO: do something
    }
    
    std::vector<uint_fast8_t> the_qualities = make_qualities(bam_get_qual(b_), sequence_length);
    
    if (the_qualities[0] == 0xff) {
        // TODO: do something
    }
    
    std::string the_sequence {make_sequence(bam_get_seq(b_), sequence_length)};
    CigarString the_cigar_string {make_cigar_string(bam_get_cigar(b_), c.n_cigar)};
    
    uint_fast32_t read_start {static_cast<uint_fast32_t>(c.pos)};
    
    // Note htslib uses 0-based positions
    return AlignedRead {
        GenomicRegion(get_contig_name(c.tid), read_start, read_start + sequence_length),
        std::move(the_sequence),
        std::move(the_qualities),
        std::move(the_cigar_string),
        static_cast<uint_fast32_t>(c.isize),
        get_contig_name(c.mtid),
        static_cast<uint_fast32_t>(c.mpos),
        static_cast<uint_fast8_t>(c.qual)
    };
}

int32_t HtslibFacade::get_htslib_tid(const std::string& reference_contig_name) const
{
    // Could avoid this lookup if stored inverse tid mappings, but I don't think this
    // function will be called much, so probably better not to explicitly store map.
    for (const auto& pair : htslib_tid_to_contig_name_) {
        if (pair.second == reference_contig_name) return pair.first;
    }
    throw std::runtime_error {"Reference contig not found"};
}
