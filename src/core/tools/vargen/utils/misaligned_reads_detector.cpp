// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "misaligned_reads_detector.hpp"

#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "basics/cigar_string.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/maths.hpp"
#include "utils/append.hpp"
#include "utils/free_memory.hpp"

namespace octopus { namespace coretools {

MisalignedReadsDetector::MisalignedReadsDetector(const ReferenceGenome& reference) : reference_ {reference} {}

MisalignedReadsDetector::MisalignedReadsDetector(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
{}

void MisalignedReadsDetector::add(const SampleName& sample, const AlignedRead& read)
{
    coverage_tracker_[sample].add(read);
    if (is_likely_misaligned(read)) {
        likely_misaligned_coverage_tracker_[sample].add(read);
    }
}

std::vector<GenomicRegion> MisalignedReadsDetector::generate(const GenomicRegion& region) const
{
    if (likely_misaligned_coverage_tracker_.empty()) return {};
    if (coverage_tracker_.size() == 1) {
        const auto& sample = std::cbegin(coverage_tracker_)->first;
        const auto total_coverages = coverage_tracker_.at(sample).get(region);
        const auto misaligned_coverages = likely_misaligned_coverage_tracker_.at(sample).get(region);
        assert(total_coverages.size() == misaligned_coverages.size());
        std::vector<bool> likely_misaligned_base_mask(total_coverages.size());
        std::transform(std::cbegin(total_coverages), std::cend(total_coverages), std::cbegin(misaligned_coverages),
                       std::begin(likely_misaligned_base_mask),
                       [] (auto depth, auto misaligned_depth) -> bool {
                           return misaligned_depth > depth / 2;
                       });
        return join(select_regions(region, likely_misaligned_base_mask), 30);
    }
    return {};
}

void MisalignedReadsDetector::clear() noexcept
{
    free_memory(coverage_tracker_);
    free_memory(likely_misaligned_coverage_tracker_);
}

namespace {

using NucleotideSequenceIterator = AlignedRead::NucleotideSequence::const_iterator;
using BaseQualityVectorIterator  = AlignedRead::BaseQualityVector::const_iterator;

bool count_snvs_in_match_range(const NucleotideSequenceIterator first_ref, const NucleotideSequenceIterator last_ref,
                               const NucleotideSequenceIterator first_base, const BaseQualityVectorIterator first_quality,
                               const AlignedRead::BaseQuality trigger)
{
    using boost::make_zip_iterator;
    using Tuple = boost::tuple<char, char, AlignedRead::BaseQuality>;
    const auto num_bases = std::distance(first_ref, last_ref);
    const auto last_base = std::next(first_base, num_bases);
    const auto last_quality = std::next(first_quality, num_bases);
    return std::count_if(make_zip_iterator(boost::make_tuple(first_ref, first_base, first_quality)),
                        make_zip_iterator(boost::make_tuple(last_ref, last_base, last_quality)),
                        [trigger](const Tuple& t) {
                            const char ref_base{t.get<0>()}, read_base{t.get<1>()};
                            return ref_base != read_base && ref_base != 'N' && read_base != 'N' && t.get<2>() >= trigger;
                        });
}

double ln_probability_read_correctly_aligned(const double misalign_penalty, const AlignedRead& read,
                                             const double max_expected_mutation_rate)
{
    const auto k = static_cast<unsigned>(std::floor(misalign_penalty));
    if (k == 0) {
        return 0;
    } else {
        const auto ln_prob_missmapped = -maths::constants::ln10Div10<> * read.mapping_quality();
        const auto ln_prob_mapped = std::log(1.0 - std::exp(ln_prob_missmapped));
        const auto mu = max_expected_mutation_rate * region_size(read);
        auto ln_prob_given_mapped = maths::log_poisson_sf(k, mu);
        return ln_prob_mapped + ln_prob_given_mapped;
    }
}

} // namespace

bool MisalignedReadsDetector::is_likely_misaligned(const AlignedRead& read) const
{
    using std::cbegin; using std::next;
    using Flag = CigarOperation::Flag;
    const auto& read_sequence = read.sequence();
    auto sequence_itr = cbegin(read_sequence);
    auto base_quality_itr = cbegin(read.base_qualities());
    auto ref_index = mapped_begin(read);
    std::size_t read_index {0};
    double misalignment_penalty {0};
    for (const auto& cigar_operation : read.cigar()) {
        const auto op_size = cigar_operation.size();
        switch (cigar_operation.flag()) {
            case Flag::alignmentMatch:
            {
                const GenomicRegion region {contig_name(read), ref_index, ref_index + op_size};
                const auto ref_segment = reference_.get().fetch_sequence(region);
                auto num_snvs = count_snvs_in_match_range(std::cbegin(ref_segment), std::cend(ref_segment),
                                                          next(sequence_itr, read_index),
                                                          next(base_quality_itr, read_index),
                                                          options_.snv_threshold);
                misalignment_penalty += num_snvs * options_.snv_penalty;
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::sequenceMatch:
                read_index += op_size;
                ref_index += op_size;
                break;
            case Flag::substitution:
            {
                auto num_snvs = std::count_if(next(base_quality_itr, read_index), next(base_quality_itr, read_index + op_size),
                                              [this] (const AlignedRead::BaseQuality quality) { return quality >= options_.snv_threshold; });
                misalignment_penalty += num_snvs * options_.snv_penalty;
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::insertion:
            {
                read_index += op_size;
                misalignment_penalty += options_.indel_penalty;
                break;
            }
            case Flag::deletion:
            {
                ref_index += op_size;
                misalignment_penalty += options_.indel_penalty;
                break;
            }
            case Flag::softClipped:
            {
                if (op_size > options_.max_unpenalised_clip_size) {
                    misalignment_penalty += options_.clip_penalty;
                }
                read_index += op_size;
                ref_index += op_size;
                break;
            }
            case Flag::hardClipped:
            {
                if (op_size > options_.max_unpenalised_clip_size) {
                    misalignment_penalty += options_.clip_penalty;
                }
                break;
            }
            case Flag::padding:
                ref_index += op_size;
                break;
            case Flag::skipped:
                ref_index += op_size;
                break;
        }
    }
    auto mu = options_.max_expected_mutation_rate;
    auto ln_prob_misaligned = ln_probability_read_correctly_aligned(misalignment_penalty, read, mu);
    auto min_ln_prob_misaligned = options_.min_ln_prob_correctly_aligned;
    return ln_prob_misaligned < min_ln_prob_misaligned;
}

} // namespace coretools
} // namespace octopus
