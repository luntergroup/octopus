// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indel_profiler_hpp
#define indel_profiler_hpp

#include <vector>
#include <deque>
#include <cstddef>
#include <iosfwd>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/aligned_read.hpp"
#include "basics/tandem_repeat.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/allele.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/variant/vcf_reader.hpp"
#include "readpipe/read_pipe.hpp"
#include "readpipe/buffered_read_pipe.hpp"
#include "utils/thread_pool.hpp"
#include "read_assigner.hpp"

namespace octopus {

class IndelProfiler
{
public:
    using SampleList = std::vector<SampleName>;
    
    struct ProfileConfig
    {
        unsigned min_period = 1, max_period = 30;
        unsigned min_periods = 1, max_periods = 50;
        unsigned max_length = 200;
        bool check_read_misalignments = true;
        Haplotype::NucleotideSequence complex_motif = "N";
    };
    
    struct PerformanceConfig
    {
        boost::optional<BufferedReadPipe::Config> read_buffering = boost::none;
        boost::optional<unsigned> max_threads = 1;
    };
    
    struct IndelProfile;
    
    IndelProfiler() = default;
    IndelProfiler(ProfileConfig config);
    IndelProfiler(ProfileConfig config, PerformanceConfig performance_config);
    
    IndelProfiler(const IndelProfiler&)            = default;
    IndelProfiler& operator=(const IndelProfiler&) = default;
    IndelProfiler(IndelProfiler&&)                 = default;
    IndelProfiler& operator=(IndelProfiler&&)      = default;
    
    ~IndelProfiler() = default;
    
    IndelProfile profile(const ReadPipe& src, VcfReader& variants, const ReferenceGenome& reference) const;
    IndelProfile profile(const ReadPipe& src, VcfReader& variants, const ReferenceGenome& reference,
                         const InputRegionMap& regions) const;
    IndelProfile profile(const ReadPipe& src, VcfReader& variants, const ReferenceGenome& reference,
                         const GenomicRegion& region) const;

private:
    using VcfIterator = VcfReader::RecordIterator;
    using CallBlock = std::vector<VcfRecord>;
    using HaplotypeReadSupportMap = std::unordered_map<Haplotype, MappableFlatMultiSet<AlignedRead>>;
    using HaplotypeRepeatMap = std::unordered_map<Haplotype, MappableFlatSet<TandemRepeat>>;
    
    struct DataBatch
    {
        Haplotype reference;
        HaplotypeReadSupportMap support;
        MappableFlatSet<Allele> indels;
        HaplotypeRepeatMap repeats;
    };
    
    ProfileConfig config_;
    PerformanceConfig performance_config_;
    mutable ThreadPool workers_;
    
    void check_samples(const SampleList& samples, const VcfReader& variants) const;
    CallBlock read_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const;
    DataBatch read_next_data_batch(VcfIterator& first, const VcfIterator& last, const ReadPipe& src,
                                   const ReferenceGenome& reference, const SampleList& samples,
                                   const GenomicRegion& analysis_region,
                                   const boost::optional<GenomicRegion>& prev_batch_region) const;
    void evaluate_indel_profile(const DataBatch& data, IndelProfile& result) const;
    void evaluate_support(const Genotype<Haplotype>& genotype, ReadContainer& reads, HaplotypeReadSupportMap& result) const;
    MappableFlatSet<TandemRepeat> find_repeats(const Haplotype& haplotype) const;
};

struct IndelProfiler::IndelProfile
{
    struct RepeatState
    {
        Haplotype::NucleotideSequence motif;
        GenomicRegion::Size span;
        unsigned reference_count, read_count;
        std::vector<unsigned> polymorphism_counts, error_counts;
    };
    using RepeaStateArray = std::vector<std::vector<std::deque<RepeatState>>>;
    RepeaStateArray states;
};

IndelProfiler::IndelProfile
profile_indels(const ReadPipe& reads, VcfReader::Path variants, const ReferenceGenome& reference);
IndelProfiler::IndelProfile
profile_indels(const ReadPipe& reads, VcfReader::Path variants, const ReferenceGenome& reference, const InputRegionMap& regions);
IndelProfiler::IndelProfile
profile_indels(const ReadPipe& reads, VcfReader::Path variants, const ReferenceGenome& reference, const GenomicRegion& region);

std::ostream& operator<<(std::ostream& os, const IndelProfiler::IndelProfile& indel_profile);

} // namespace octopus

#endif
