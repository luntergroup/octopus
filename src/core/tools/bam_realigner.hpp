// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef bam_realigner_hpp
#define bam_realigner_hpp

#include <vector>
#include <cstddef>

#include <boost/optional.hpp>

#include "basics/aligned_read.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "containers/mappable_flat_set.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/read/read_reader.hpp"
#include "io/read/read_writer.hpp"
#include "io/variant/vcf_reader.hpp"
#include "utils/memory_footprint.hpp"
#include "utils/thread_pool.hpp"

namespace octopus {

class BAMRealigner
{
public:
    using ReadReader = io::ReadReader;
    using ReadWriter = io::ReadWriter;
    using SampleName = ReadReader::SampleName;
    using SampleList = std::vector<SampleName>;
    
    struct Config
    {
        HaplotypeLikelihoodModel alignment_model = {};
        bool primary_only = true;
        bool copy_hom_ref_reads = false;
        bool simplify_cigars = false;
        bool use_paired_reads = false;
        bool use_linked_reads = false;
        MemoryFootprint max_buffer = *parse_footprint("50M");
        boost::optional<unsigned> max_threads = 1;
    };
    
    struct Report
    {
        std::size_t n_reads_assigned;
        std::size_t n_reads_unassigned;
    };
    
    BAMRealigner() = default;
    BAMRealigner(Config config);
    
    BAMRealigner(const BAMRealigner&)            = delete;
    BAMRealigner& operator=(const BAMRealigner&) = delete;
    BAMRealigner(BAMRealigner&&)                 = delete;
    BAMRealigner& operator=(BAMRealigner&&)      = delete;
    
    ~BAMRealigner() = default;
    
    Report realign(ReadReader& src, VcfReader& variants, ReadWriter& dst,
                   const ReferenceGenome& reference, SampleList samples) const;
    Report realign(ReadReader& src, VcfReader& variants, ReadWriter& dst,
                   const ReferenceGenome& reference) const;
    
private:
    using VcfIterator = VcfReader::RecordIterator;
    using CallBlock   = std::vector<VcfRecord>;
    struct Batch
    {
        MappableFlatSet<Genotype<Haplotype>> genotypes;
        std::vector<AlignedRead> reads;
    };
    using BatchList = std::vector<Batch>;
    using BatchListRegionPair = std::pair<BatchList, boost::optional<GenomicRegion>>;
    
    Config config_;
    mutable ThreadPool workers_;
    
    CallBlock read_next_block(VcfIterator& first, const VcfIterator& last, const SampleList& samples) const;
    BatchListRegionPair read_next_batch(VcfIterator& first, const VcfIterator& last, ReadReader& src,
                                        const ReferenceGenome& reference, const SampleList& samples,
                                        const boost::optional<GenomicRegion>& prev_batch_region) const;
    void merge(BatchList& src, BatchList& dst) const;
};

BAMRealigner::Report
realign(io::ReadReader::Path src, VcfReader::Path variants, io::ReadWriter::Path dst,
        const ReferenceGenome& reference);
BAMRealigner::Report
realign(io::ReadReader::Path src, VcfReader::Path variants, io::ReadWriter::Path dst,
        const ReferenceGenome& reference, BAMRealigner::Config config);

} // namespace octopus

#endif
