// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "dynamic_cigar_scanner.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "config/common.hpp"
#include "basics/aligned_read.hpp"
#include "basics/cigar_string.hpp"
#include "io/reference/reference_genome.hpp"
#include "concepts/mappable_range.hpp"
#include "utils/mappable_algorithms.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace coretools {

std::unique_ptr<VariantGenerator> DynamicCigarScanner::do_clone() const
{
    return std::make_unique<DynamicCigarScanner>(*this);
}

DynamicCigarScanner::DynamicCigarScanner(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
, candidates_ {}
, max_seen_candidate_size_ {}
, read_coverage_tracker_ {}
{}

bool DynamicCigarScanner::do_requires_reads() const noexcept
{
    return true;
}

namespace {

template <typename Sequence, typename P, typename S>
Sequence splice(const Sequence& sequence, const P pos, const S size)
{
    const auto it = std::next(std::cbegin(sequence), pos);
    return Sequence {it, std::next(it, size)};
}

} // namespace

void DynamicCigarScanner::do_add_read(const AlignedRead& read)
{
    using std::cbegin; using std::next; using std::move;
    using Flag = CigarOperation::Flag;
    
    const auto& read_contig   = contig_name(read);
    const auto& read_sequence = read.sequence();
    auto sequence_iter     = cbegin(read_sequence);
    auto base_quality_iter = cbegin(read.qualities());
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
    GenomicRegion region;
    
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        switch (cigar_operation.flag()) {
        case Flag::alignmentMatch:
            add_snvs_in_match_range(GenomicRegion {read_contig, ref_index, ref_index + op_size},
                                    next(sequence_iter, read_index),
                                    next(sequence_iter, read_index + op_size),
                                    next(base_quality_iter, read_index));
            read_index += op_size;
            ref_index  += op_size;
            break;
        case Flag::sequenceMatch:
            read_index += op_size;
            ref_index  += op_size;
            break;
        case Flag::substitution:
        {
            region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
            add_candidate(region,
                          reference_.get().fetch_sequence(region),
                          splice(read_sequence, read_index, op_size),
                          next(base_quality_iter, read_index));
            read_index += op_size;
            ref_index  += op_size;
            break;
        }
        case Flag::insertion:
        {
            add_candidate(GenomicRegion {read_contig, ref_index, ref_index},
                          "",
                          splice(read_sequence, read_index, op_size),
                          next(base_quality_iter, read_index));
            read_index += op_size;
            break;
        }
        case Flag::deletion:
        {
            region = GenomicRegion {read_contig, ref_index, ref_index + op_size};
            add_candidate(move(region),
                          reference_.get().fetch_sequence(region),
                          "",
                          next(base_quality_iter, read_index));
            ref_index += op_size;
            break;
        }
        case Flag::softClipped:
            read_index += op_size;
            ref_index  += op_size;
            break;
        case Flag::hardClipped:
            break;
        case Flag::padding:
            ref_index += op_size;
            break;
        case Flag::skipped:
            ref_index += op_size;
            break;
        }
    }
    if (options_.use_clipped_coverage_tracking) {
        read_coverage_tracker_.add(clipped_mapped_region(read));
    } else {
        read_coverage_tracker_.add(read);
    }
}

void DynamicCigarScanner::do_add_reads(VectorIterator first, VectorIterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { do_add_read(read); });
}

void DynamicCigarScanner::do_add_reads(FlatSetIterator first, FlatSetIterator last)
{
    std::for_each(first, last, [this] (const auto& read ) { do_add_read(read); });
}

std::vector<Variant> DynamicCigarScanner::do_generate_variants(const GenomicRegion& region)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::next; using std::distance;
    
    std::sort(begin(candidates_), end(candidates_));
    auto overlapped = overlap_range(candidates_, region, max_seen_candidate_size_);
    
    std::vector<Variant> result {};
    result.reserve(size(overlapped, BidirectionallySortedTag {})); // maximum possible
    
    while (!overlapped.empty()) {
        const Candidate& candidate {overlapped.front()};
        const auto next_candidate = std::find_if_not(next(cbegin(overlapped)), cend(overlapped),
                                                     [this, &candidate] (const Candidate& c) {
                                                         return options_.match(c.variant, candidate.variant);
                                                     });
        const auto num_observations = static_cast<unsigned>(distance(cbegin(overlapped), next_candidate));
        const auto min_depth = read_coverage_tracker_.min_coverage(mapped_region(candidate));
        const auto base_quality_sum = std::accumulate(cbegin(overlapped), next_candidate, 0u,
                                                      [this] (auto curr, const Candidate& candidate) {
                                                          return curr + sum_base_qualities(candidate);
                                                      });
        if (options_.include(candidate.variant, num_observations, min_depth, base_quality_sum)) {
            if (num_observations > 1) {
                auto unique_iter = cbegin(overlapped);
                while (unique_iter != next_candidate) {
                    result.push_back(unique_iter->variant);
                    unique_iter = std::find_if_not(next(unique_iter), next_candidate,
                                                   [unique_iter] (const Candidate& c) {
                                                       return c.variant == unique_iter->variant;
                                                   });
                }
            } else {
                result.push_back(candidate.variant);
            }
        }
        overlapped.advance_begin(num_observations);
    }
    
    return result;
}

void DynamicCigarScanner::do_clear() noexcept
{
    candidates_.clear();
    candidates_.shrink_to_fit();
}

std::string DynamicCigarScanner::name() const
{
    return "DynamicCigarScanner";
}

// private methods

void DynamicCigarScanner::add_snvs_in_match_range(const GenomicRegion& region,
                                                  const SequenceIterator first_base, const SequenceIterator last_base,
                                                  AlignedRead::BaseQualityVector::const_iterator first_base_quality)
{
    using boost::make_zip_iterator; using std::for_each; using std::cbegin; using std::cend;
    using Tuple = boost::tuple<char, char>;
    
    const NucleotideSequence ref_segment {reference_.get().fetch_sequence(region)};
    const auto& contig = region.contig_name();
    auto ref_index = mapped_begin(region);
    
    for_each(make_zip_iterator(boost::make_tuple(cbegin(ref_segment), first_base)),
             make_zip_iterator(boost::make_tuple(cend(ref_segment), last_base)),
             [this, &contig, &ref_index, &first_base_quality] (const Tuple& t) {
                 const char ref_base  {t.get<0>()}, read_base {t.get<1>()};
                 if (ref_base != read_base && ref_base != 'N' && read_base != 'N') {
                     add_candidate(GenomicRegion {contig, ref_index, ref_index + 1},
                                   ref_base, read_base, first_base_quality);
                 }
                 ++ref_index;
                 ++first_base_quality;
             });
}

unsigned DynamicCigarScanner::sum_base_qualities(const Candidate& candidate) const noexcept
{
    return std::accumulate(candidate.first_base_quality_iter,
                           std::next(candidate.first_base_quality_iter, alt_sequence_size(candidate.variant)),
                           0u);
}

} // coretools
} // namespace octopus
