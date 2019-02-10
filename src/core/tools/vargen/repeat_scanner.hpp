// Copyright (c) 2015-2019 Daniel Cooke
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
        unsigned min_observations = 2;
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
    std::vector<Variant> do_generate(const RegionSet& regions) const override;
    void do_clear() noexcept override;
    std::string name() const override;
    
    struct SNV : public Comparable<SNV>, public Mappable<SNV>
    {
        SNV(ContigRegion::Position position, char base) : region {position, position + 1}, base {base} {}
        ContigRegion region;
        char base;
        const ContigRegion& mapped_region() const noexcept { return region; }
    };
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Options options_;
    
    mutable std::vector<SNV> snv_buffer_;
    mutable std::deque<Variant> candidates_;
    
    void add_read(const AlignedRead& read);
    void add_match_range(const GenomicRegion& region, const AlignedRead& read, std::size_t read_index) const;
    void add_to_buffer(SNV snv, const ContigName& contig) const;
    void reset_buffer(const ContigName& contig) const;
    void generate(const GenomicRegion& region, std::vector<Variant>& result) const;
    std::vector<Variant> get_candidate_mnvs(const GenomicRegion& region) const;
};

} // coretools
} // namespace octopus

#endif
