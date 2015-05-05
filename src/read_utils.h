//
//  read_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_utils__
#define __Octopus__read_utils__

#include <vector>
#include <unordered_map>
#include <iterator> // std::begin, std::end, std::make_move_iterator
#include <utility>  // std::move

#include "aligned_read.h"
#include "read_filter.h"

template <typename T, typename Container>
using ReadMap = std::unordered_map<T, Container>;

template <typename Container>
void reserve_if_enabled(Container& container, typename Container::size_type n) {}

template <typename T>
void reserve_if_enabled(std::vector<T>& container, typename std::vector<T>::size_type n) { container.reserve(n); }

template <typename Container>
void shrink_to_fit_if_enabled(Container& container) {}

template <typename T>
void shrink_to_fit_if_enabled(std::vector<T>& container) { container.shrink_to_fit(); }

template <typename T, typename Container, typename ReadFilter>
std::pair<ReadMap<T, Container>, ReadMap<T, Container>>
filter_reads(ReadMap<T, Container>&& the_reads, ReadFilter& a_read_filter)
{
    ReadMap<T, Container> good_read_map {}, bad_read_map {};
    good_read_map.reserve(the_reads.size());
    bad_read_map.reserve(the_reads.size());
    
    for (auto& sample_reads : the_reads) {
        Container good_reads {}, bad_reads {};
        
        reserve_if_enabled(good_reads, sample_reads.second.size());
        reserve_if_enabled(bad_reads, sample_reads.second.size() / 10); // arbitrarily chosen
        
        a_read_filter.filter_reads(std::make_move_iterator(std::begin(sample_reads.second)),
                                   std::make_move_iterator(std::end(sample_reads.second)),
                                   ContextBackInserter(good_reads), ContextBackInserter(bad_reads));
        
        sample_reads.second.clear();
        shrink_to_fit_if_enabled(good_reads);
        shrink_to_fit_if_enabled(bad_reads);
        
        good_read_map.emplace(sample_reads.first, std::move(good_reads));
        bad_read_map.emplace(std::move(sample_reads.first), std::move(bad_reads));
    }
    
    return {good_read_map, bad_read_map};
}

template <typename ForwardIterator>
unsigned get_min_coverage(ForwardIterator first, ForwardIterator last)
{
    return 0;
}

unsigned get_min_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

unsigned get_max_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

double get_mean_coverage(const std::vector<AlignedRead>& reads, const GenomicRegion& a_region);

#endif /* defined(__Octopus__read_utils__) */
