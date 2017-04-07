// Copyright (c) 2016 Daniel Cooke
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
    struct Options
    {
        // Signature is (Variant, depth, observed base_qualities). For insertions, observed base_qualities is the sum.
        using InclusionPredicate = std::function<bool(const Variant&, unsigned, std::vector<unsigned>&)>;
        using MatchPredicate = std::function<bool(const Variant&, const Variant&)>;
        using RepeatRegionGenerator = std::function<std::vector<GenomicRegion>(const ReferenceGenome&, GenomicRegion)>;
        InclusionPredicate include;
        MatchPredicate match = std::equal_to<> {};
        bool use_clipped_coverage_tracking = false;
        Variant::MappingDomain::Size max_variant_size = 2000;
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
    void add_read(const AlignedRead& read);
    void do_add_read(const SampleName& sample, const AlignedRead& read) override;
    void do_add_reads(const SampleName& sample, VectorIterator first, VectorIterator last) override;
    void do_add_reads(const SampleName& sample, FlatSetIterator first, FlatSetIterator last) override;
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    void do_clear() noexcept override;
    std::string name() const override;
    
    struct Candidate : public Comparable<Candidate>, public Mappable<Candidate>
    {
        Variant variant;
        AlignedRead::BaseQualityVector::const_iterator first_base_quality_iter;
        template <typename T1, typename T2, typename T3>
        Candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                  AlignedRead::BaseQualityVector::const_iterator first_base_quality);
        const GenomicRegion& mapped_region() const noexcept { return variant.mapped_region(); }
        friend bool operator==(const Candidate& lhs, const Candidate& rhs) noexcept { return lhs.variant == rhs.variant; }
        friend bool operator<(const Candidate& lhs, const Candidate& rhs) noexcept { return lhs.variant < rhs.variant; }
    };
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using SequenceIterator   = NucleotideSequence::const_iterator;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Options options_;
    std::deque<Candidate> candidates_;
    Variant::MappingDomain::Size max_seen_candidate_size_;
    CoverageTracker<GenomicRegion> read_coverage_tracker_;
    
    template <typename T1, typename T2, typename T3>
    void add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                       AlignedRead::BaseQualityVector::const_iterator first_base_quality);
    void add_snvs_in_match_range(const GenomicRegion& region,
                                 SequenceIterator first_base, SequenceIterator last_base,
                                 AlignedRead::BaseQualityVector::const_iterator first_quality);
    unsigned sum_base_qualities(const Candidate& candidate) const noexcept;
    std::vector<GenomicRegion> get_repeat_regions(const GenomicRegion& region) const;
};

template <typename T1, typename T2, typename T3>
DynamicCigarScanner::Candidate::Candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                                          AlignedRead::BaseQualityVector::const_iterator first_base_quality)
: variant {std::forward<T1>(region), std::forward<T2>(sequence_removed), std::forward<T3>(sequence_added)}
, first_base_quality_iter {first_base_quality}
{}

template <typename T1, typename T2, typename T3>
void DynamicCigarScanner::add_candidate(T1&& region, T2&& sequence_removed, T3&& sequence_added,
                                        AlignedRead::BaseQualityVector::const_iterator first_base_quality)
{
    const auto candidate_size = size(region);
    if (candidate_size <= options_.max_variant_size) {
        candidates_.emplace_back(std::forward<T1>(region),
                                 std::forward<T2>(sequence_removed),
                                 std::forward<T3>(sequence_added),
                                 first_base_quality);
        max_seen_candidate_size_ = std::max(max_seen_candidate_size_, candidate_size);
    }
}

} // coretools
} // namespace octopus

#endif
