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
    for (const auto& pair : htslib_tid_to_contig_name_) {
        if (pair.second == reference_contig_name) return pair.first;
    }
    throw std::runtime_error {"Reference contig not found"};
}
