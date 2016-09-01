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
        Variant::MappingDomain::Size max_variant_size = 100;
        bool always_include_overlapping_indels        = true;
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
    
    struct Candidate
    {
        Variant variant;
        AlignedRead::BaseQuality quality;
    };
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    
    bool do_requires_reads() const noexcept override;
    
    void do_add_read(const AlignedRead& read) override;
    void do_add_reads(VectorIterator first, VectorIterator last) override;
    void do_add_reads(FlatSetIterator first, FlatSetIterator last) override;
    
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    
    void do_clear() noexcept override;
    
    std::string name() const override;
    
    using NucleotideSequence = AlignedRead::NucleotideSequence;
    using SequenceIterator   = NucleotideSequence::const_iterator;
    using QualitiesIterator  = AlignedRead::BaseQualityVector::const_iterator;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    Options options_;
    
    std::function<bool(const Variant&, const Variant&)> match_;
    
    std::deque<Variant> candidates_;
    
    Variant::MappingDomain::Size max_seen_candidate_size_;
    
    CoverageTracker read_coverage_tracker_;
    
    void add_snvs_in_match_range(const GenomicRegion& region, SequenceIterator first_base,
                                 SequenceIterator last_base, QualitiesIterator first_quality);
};

} // coretools
} // namespace octopus

#endif
