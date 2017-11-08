// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef dynamic_cigar_scanner_hpp
#define dynamic_cigar_scanner_hpp

#include <vector>
#include <deque>
#include <cstddef>
#include <utility>
#include <functional>
#include <memory>

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

class DynamicCigarScanner : public VariantGenerator
{
public:
    struct ObservedVariant
    {
        Variant variant;
        unsigned num_samples, total_depth;
        struct SampleObservation
        {
            unsigned depth;
            std::vector<unsigned> observed_qualities;
            unsigned num_fwd_observations;
        };
        std::vector<SampleObservation> sample_observations;
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
        
        using InclusionPredicate = std::function<bool(ObservedVariant)>;
        using MatchPredicate = std::function<bool(const Variant&, const Variant&)>;
        using RepeatRegionGenerator = std::function<std::vector<GenomicRegion>(const ReferenceGenome&, GenomicRegion)>;
        InclusionPredicate include;
        MatchPredicate match = std::equal_to<> {};
        bool use_clipped_coverage_tracking = false;
        Variant::MappingDomain::Size max_variant_size = 2000;
        MisalignmentParameters misalignment_parameters = MisalignmentParameters {};
        boost::optional<RepeatRegionGenerator> repeat_region_generator = boost::none;
        double max_repeat_region_density = 2;
    };
    
    DynamicCigarScanner() = delete;
    DynamicCigarScanner(const ReferenceGenome& reference, Options options);
    
    DynamicCigarScanner(const DynamicCigarScanner&)            = default;
    DynamicCigarScanner& operator=(const DynamicCigarScanner&) = default;
    DynamicCigarScanner(DynamicCigarScanner&&)                 = default;
    DynamicCigarScanner& operator=(DynamicCigarScanner&&)      = default;
    
    ~DynamicCigarScanner() override = default;
    
private:
    using VariantGenerator::VectorIterator;
    using VariantGenerator::FlatSetIterator;
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    bool do_requires_reads() const noexcept override;
    void do_add_read(const SampleName& sample, const AlignedRead& read) override;
    void add_read(const SampleName& sample, const AlignedRead& read,
                  CoverageTracker<GenomicRegion>& sample_coverage_tracker);
    void do_add_reads(const SampleName& sample, VectorIterator first, VectorIterator last) override;
    void do_add_reads(const SampleName& sample, FlatSetIterator first, FlatSetIterator last) override;
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    void do_clear() noexcept override;
    std::string name() const override;
    
    struct Candidate : public Comparable<Candidate>, public Mappable<Candidate>
    {
        Variant variant;
        SampleName origin;
        AlignedRead::BaseQualityVector::const_iterator first_base_quality_iter;
        AlignedRead::Direction support_direction;
        
        template <typename T1, typename T2, typename T3, typename T4>
        Candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added, T4&& origin,
                  AlignedRead::BaseQualityVector::const_iterator first_base_quality,
                  AlignedRead::Direction support_direction);
        
        const GenomicRegion& mapped_region() const noexcept { return variant.mapped_region(); }
        
        friend bool operator==(const Candidate& lhs, const Candidate& rhs) noexcept { return lhs.variant == rhs.variant; }
        friend bool operator<(const Candidate& lhs, const Candidate& rhs) noexcept { return lhs.variant < rhs.variant; }
    };
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using SequenceIterator   = NucleotideSequence::const_iterator;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Options options_;
    std::vector<Candidate> buffer_;
    std::deque<Candidate> candidates_, likely_misaligned_candidates_;
    Variant::MappingDomain::Size max_seen_candidate_size_;
    CoverageTracker<GenomicRegion> read_coverage_tracker_, misaligned_tracker_;
    std::unordered_map<SampleName, CoverageTracker<GenomicRegion>> sample_read_coverage_tracker_;
    
    template <typename T1, typename T2, typename T3, typename T4>
    void add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                       T4&& origin,
                       AlignedRead::BaseQualityVector::const_iterator first_base_quality,
                       AlignedRead::Direction support_direction);
    double add_snvs_in_match_range(const GenomicRegion& region, SequenceIterator first_base, SequenceIterator last_base,
                                   const SampleName& origin,
                                   AlignedRead::BaseQualityVector::const_iterator first_quality,
                                   AlignedRead::Direction support_direction);
    unsigned sum_base_qualities(const Candidate& candidate) const noexcept;
    std::vector<GenomicRegion> get_repeat_regions(const GenomicRegion& region) const;
};

template <typename T1, typename T2, typename T3, typename T4>
DynamicCigarScanner::Candidate::Candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                                          T4&& origin,
                                          AlignedRead::BaseQualityVector::const_iterator first_base_quality,
                                          AlignedRead::Direction support_direction)
: variant {std::forward<T1>(region), std::forward<T2>(sequence_removed), std::forward<T3>(sequence_added)}
, origin {std::forward<T4>(origin)}
, first_base_quality_iter {first_base_quality}
, support_direction {support_direction}
{}

template <typename T1, typename T2, typename T3, typename T4>
void DynamicCigarScanner::add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                                        T4&& origin,
                                        AlignedRead::BaseQualityVector::const_iterator first_base_quality,
                                        AlignedRead::Direction support_direction)
{
    const auto candidate_size = size(region);
    if (candidate_size <= options_.max_variant_size) {
        buffer_.emplace_back(std::forward<T1>(region),
                             std::forward<T2>(sequence_removed),
                             std::forward<T3>(sequence_added),
                             std::forward<T4>(origin),
                             first_base_quality,
                             support_direction);
        max_seen_candidate_size_ = std::max(max_seen_candidate_size_, candidate_size);
    }
}

struct DefaultInclusionPredicate
{
    bool operator()(DynamicCigarScanner::ObservedVariant);
};

struct DefaultMatchPredicate
{
    bool operator()(const Variant& lhs, const Variant& rhs) noexcept;
};

} // coretools
} // namespace octopus

#endif
