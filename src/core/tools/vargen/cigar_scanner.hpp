// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cigar_scanner_hpp
#define cigar_scanner_hpp

#include <vector>
#include <deque>
#include <cstddef>
#include <utility>
#include <functional>
#include <memory>

#include <boost/optional.hpp>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/variant.hpp"
#include "utils/coverage_tracker.hpp"
#include "variant_generator.hpp"

namespace octopus {

class ReferenceGenome;
class GenomicRegion;

namespace coretools {

class CigarScanner : public VariantGenerator
{
public:
    struct VariantObservation
    {
        Variant variant;
        unsigned total_depth;
        struct SampleObservationStats
        {
            std::reference_wrapper<const SampleName> sample;
            unsigned depth, forward_strand_depth;
            std::vector<unsigned> observed_base_qualities;
            std::vector<AlignedRead::MappingQuality> observed_mapping_qualities;
            unsigned forward_strand_support, edge_support;
        };
        std::vector<SampleObservationStats> sample_observations;
    };
    
    struct Options
    {
        struct MisalignmentParameters
        {
            AlignedRead::BaseQuality snv_threshold;
            double snv_penalty = 1, indel_penalty = 1, clip_penalty = 1;
            double max_expected_mutation_rate = 1e-3;
            double min_ln_prob_correctly_aligned = std::log(0.0001);
            unsigned max_unpenalised_clip_size = 3;
        };
        
        using InclusionPredicate = std::function<bool(VariantObservation)>;
        using MatchPredicate = std::function<bool(const Variant&, const Variant&)>;
        
        InclusionPredicate include;
        MatchPredicate match = std::equal_to<> {};
        bool use_clipped_coverage_tracking = false;
        Variant::MappingDomain::Size max_variant_size = 2000;
        boost::optional<MisalignmentParameters> misalignment_parameters = MisalignmentParameters {};
        bool split_mnvs = true;
        bool ignore_strand_bias = false;
    };
    
    CigarScanner() = delete;
    CigarScanner(const ReferenceGenome& reference, Options options);
    
    CigarScanner(const CigarScanner&)            = default;
    CigarScanner& operator=(const CigarScanner&) = default;
    CigarScanner(CigarScanner&&)                 = default;
    CigarScanner& operator=(CigarScanner&&)      = default;
    
    ~CigarScanner() override = default;
    
private:
    using VariantGenerator::ReadVectorIterator;
    using VariantGenerator::ReadFlatSetIterator;
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    bool do_requires_reads() const noexcept override;
    void do_add_read(const SampleName& sample, const AlignedRead& read) override;
    void do_add_template(const SampleName& sample, const AlignedTemplate& reads) override;
    void add_read(const SampleName& sample, const AlignedRead& read,
                  CoverageTracker<GenomicRegion>& coverage_tracker,
                  CoverageTracker<GenomicRegion>& forward_strand_coverage_tracker);
    void add_template(const SampleName& sample, const AlignedTemplate& reads,
                      CoverageTracker<GenomicRegion>& coverage_tracker,
                      CoverageTracker<GenomicRegion>& forward_strand_coverage_tracker);
    void do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last) override;
    void do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last) override;
    void do_add_reads(const SampleName& sample, TemplateVectorIterator first, TemplateVectorIterator last) override;
    void do_add_reads(const SampleName& sample, TemplateFlatSetIterator first, TemplateFlatSetIterator last) override;
    std::vector<Variant> do_generate(const RegionSet& regions, OptionalThreadPool workers) const override;
    void do_clear() noexcept override;
    std::string name() const override;
    
    struct Candidate : public Comparable<Candidate>, public Mappable<Candidate>
    {
        Variant variant;
        std::reference_wrapper<const AlignedRead> source;
        std::reference_wrapper<const SampleName> origin;
        std::size_t offset;
        
        template <typename T1, typename T2, typename T3>
        Candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                  const AlignedRead& source, std::size_t offset,
                  const SampleName& origin);
        
        const GenomicRegion& mapped_region() const noexcept { return variant.mapped_region(); }
        
        friend bool operator==(const Candidate& lhs, const Candidate& rhs) noexcept { return lhs.variant == rhs.variant; }
        friend bool operator<(const Candidate& lhs, const Candidate& rhs) noexcept { return lhs.variant < rhs.variant; }
    };
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using SequenceIterator = NucleotideSequence::const_iterator;
    using SampleCoverageTrackerMap = std::unordered_map<SampleName, CoverageTracker<GenomicRegion>>;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Options options_;
    std::vector<Candidate> buffer_;
    mutable std::deque<Candidate> candidates_, likely_misaligned_candidates_;
    Variant::MappingDomain::Size max_seen_candidate_size_;
    CoverageTracker<GenomicRegion> combined_read_coverage_tracker_, misaligned_read_coverage_tracker_;
    SampleCoverageTrackerMap sample_read_coverage_tracker_, sample_forward_strand_coverage_tracker_;
    std::deque<AlignedRead> artificial_read_buffer_;
    
    using CandidateIterator = OverlapIterator<decltype(candidates_)::const_iterator>;
    
    template <typename T1, typename T2, typename T3>
    void add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                       const AlignedRead& read, std::size_t offset, const SampleName& sample);
    double add_snvs_in_match_range(const GenomicRegion& region, const AlignedRead& read,
                                   std::size_t read_index, const SampleName& origin);
    void generate(const GenomicRegion& region, std::vector<Variant>& result) const;
    unsigned sum_base_qualities(const Candidate& candidate) const noexcept;
    bool is_likely_misaligned(const AlignedRead& read, double penalty) const;
    VariantObservation make_observation(CandidateIterator first_match, CandidateIterator last_match) const;
    std::vector<Variant> get_novel_likely_misaligned_candidates(const std::vector<Variant>& current_candidates) const;
};

template <typename T1, typename T2, typename T3>
CigarScanner::Candidate::Candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                                   const AlignedRead& source, std::size_t offset,
                                   const SampleName& origin)
: variant {std::forward<T1>(region), std::forward<T2>(sequence_removed), std::forward<T3>(sequence_added)}
, source {source}
, origin {origin}
, offset {offset}
{}

template <typename T1, typename T2, typename T3>
void CigarScanner::add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                                 const AlignedRead& read, const std::size_t offset,
                                 const SampleName& sample)
{
    const auto candidate_size = size(region);
    if (candidate_size <= options_.max_variant_size) {
        buffer_.emplace_back(std::forward<T1>(region),
                             std::forward<T2>(sequence_removed),
                             std::forward<T3>(sequence_added),
                             read, offset, sample);
        max_seen_candidate_size_ = std::max(max_seen_candidate_size_, candidate_size);
    }
}

struct KnownCopyNumberInclusionPredicate
{
    KnownCopyNumberInclusionPredicate(unsigned copy_number = 2) : copy_number_ {copy_number} {}
    bool operator()(const CigarScanner::VariantObservation& candidate);
private:
    unsigned copy_number_;
};

struct PacBioInclusionPredicate
{
    bool operator()(const CigarScanner::VariantObservation& candidate);
};

struct UnknownCopyNumberInclusionPredicate
{
    UnknownCopyNumberInclusionPredicate() = default;
    
    UnknownCopyNumberInclusionPredicate(double min_vaf, double min_probability = 0.5);
    UnknownCopyNumberInclusionPredicate(SampleName normal, double min_vaf, double min_probability = 0.5);
    
    bool operator()(const CigarScanner::VariantObservation& candidate);
private:
    boost::optional<SampleName> normal_;
    double min_vaf_ = 0.01, min_probability_ = 0.5;
};

struct CellInclusionPredicate
{
    bool operator()(const CigarScanner::VariantObservation& candidate);
};

struct SimpleThresholdInclusionPredicate
{
    SimpleThresholdInclusionPredicate(std::size_t min_observations) : min_observations_ {min_observations} {}
    bool operator()(const CigarScanner::VariantObservation& candidate) noexcept;
private:
    std::size_t min_observations_;
};

struct TolerantMatchPredicate
{
    bool operator()(const Variant& lhs, const Variant& rhs) noexcept;
};

} // coretools
} // namespace octopus

#endif
