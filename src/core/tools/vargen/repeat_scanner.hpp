// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef repeat_scanner_hpp
#define repeat_scanner_hpp

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

class RepeatScanner : public VariantGenerator
{
public:
    struct Options
    {
        unsigned min_snvs = 2;
        unsigned max_period = 6;
        unsigned min_tract_length = 3;
        unsigned min_observations = 2;
        unsigned min_sample_observations = 2;
        boost::optional<double> min_vaf = boost::none;
        AlignedRead::BaseQuality min_base_quality = 10;
    };
    
    RepeatScanner() = delete;
    RepeatScanner(const ReferenceGenome& reference, Options options);
    
    RepeatScanner(const RepeatScanner&)            = default;
    RepeatScanner& operator=(const RepeatScanner&) = default;
    RepeatScanner(RepeatScanner&&)                 = default;
    RepeatScanner& operator=(RepeatScanner&&)      = default;
    
    ~RepeatScanner() override = default;

private:
    using VariantGenerator::ReadVectorIterator;
    using VariantGenerator::ReadFlatSetIterator;
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    bool do_requires_reads() const noexcept override;
    void do_add_read(const SampleName& sample, const AlignedRead& read) override;
    void do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last) override;
    void do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last) override;
    std::vector<Variant> do_generate(const RegionSet& regions, OptionalThreadPool workers) const override;
    void do_clear() noexcept override;
    std::string name() const override;
    
    struct SNV : public Comparable<SNV>, public Mappable<SNV>
    {
        SNV(ContigRegion::Position position, char base) : region {position, position + 1}, base {base} {}
        ContigRegion region;
        char base;
        const ContigRegion& mapped_region() const noexcept { return region; }
    };
    
    struct Candidate : public Comparable<Candidate>, public Mappable<Candidate>
    {
        Candidate(Variant variant, unsigned sample_index) : variant {std::move(variant)}, sample_index {sample_index} {}
        Variant variant;
        unsigned sample_index;
        const GenomicRegion& mapped_region() const noexcept { return variant.mapped_region(); }
        friend bool operator<(const Candidate& lhs, const Candidate& rhs) { return lhs.variant < rhs.variant; }
    };
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Options options_;
    std::vector<SampleName> samples_;
    std::unordered_map<SampleName, CoverageTracker<GenomicRegion>> read_coverage_tracker_;
    
    mutable std::vector<SNV> snv_buffer_;
    mutable std::deque<Candidate> candidates_;
    
    unsigned get_sample_index(const SampleName& sample);
    void add_read(const AlignedRead& read, unsigned sample_index);
    void add_match_range(const GenomicRegion& region, const AlignedRead& read, std::size_t read_index, unsigned sample_index) const;
    void add_to_buffer(SNV snv, const ContigName& contig, unsigned sample_index) const;
    void reset_buffer(const ContigName& contig, unsigned sample_index) const;
    void generate(const GenomicRegion& region, std::vector<Variant>& result) const;
    std::vector<Variant> get_candidate_mnvs(const GenomicRegion& region) const;
};

} // coretools
} // namespace octopus

#endif
