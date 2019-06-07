// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef calling_components_hpp
#define calling_components_hpp

#include <vector>
#include <cstddef>
#include <functional>
#include <memory>

#include <boost/optional.hpp>
#include <boost/filesystem/path.hpp>

#include "config/common.hpp"
#include "config/option_parser.hpp"
#include "basics/genomic_region.hpp"
#include "basics/ploidy_map.hpp"
#include "basics/pedigree.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/read/read_manager.hpp"
#include "io/variant/vcf_writer.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "core/callers/caller_factory.hpp"
#include "core/csr/filters/variant_call_filter_factory.hpp"
#include "core/tools/bam_realigner.hpp"
#include "utils/memory_footprint.hpp"
#include "utils/input_reads_profiler.hpp"
#include "logging/progress_meter.hpp"

namespace octopus {

class GenomeCallingComponents
{
public:
    using Path = boost::filesystem::path;
    
    GenomeCallingComponents() = delete;
    
    GenomeCallingComponents(ReferenceGenome&& reference, ReadManager&& read_manager,
                            VcfWriter&& output, const options::OptionMap& options);
    
    GenomeCallingComponents(const GenomeCallingComponents&)            = delete;
    GenomeCallingComponents& operator=(const GenomeCallingComponents&) = delete;
    GenomeCallingComponents(GenomeCallingComponents&& other) noexcept;
    GenomeCallingComponents& operator=(GenomeCallingComponents&& other) = delete;
    
    ~GenomeCallingComponents() = default;
    
    const ReferenceGenome& reference() const noexcept;
    ReadManager& read_manager() noexcept;
    const ReadManager& read_manager() const noexcept;
    ReadPipe& read_pipe() noexcept;
    const ReadPipe& read_pipe() const noexcept;
    const std::vector<SampleName>& samples() const noexcept;
    const InputRegionMap& search_regions() const noexcept;
    const std::vector<GenomicRegion::ContigName>& contigs() const noexcept;
    VcfWriter& output() noexcept;
    const VcfWriter& output() const noexcept;
    MemoryFootprint read_buffer_footprint() const noexcept;
    std::size_t read_buffer_size() const noexcept;
    const boost::optional<Path>& temp_directory() const noexcept;
    boost::optional<unsigned> num_threads() const noexcept;
    const HaplotypeLikelihoodModel& haplotype_likelihood_model() const noexcept;
    const CallerFactory& caller_factory() const noexcept;
    boost::optional<VcfWriter&> filtered_output() noexcept;
    boost::optional<const VcfWriter&> filtered_output() const noexcept;
    const VariantCallFilterFactory& call_filter_factory() const;
    ReadPipe& filter_read_pipe() noexcept;
    const ReadPipe& filter_read_pipe() const noexcept;
    ProgressMeter& progress_meter() noexcept;
    bool sites_only() const noexcept;
    const PloidyMap& ploidies() const noexcept;
    boost::optional<Pedigree> pedigree() const;
    boost::optional<Path> legacy() const;
    boost::optional<Path> filter_request() const;
    boost::optional<Path> bamout() const;
    BAMRealigner::Config bamout_config() const noexcept;
    boost::optional<ReadSetProfile> reads_profile() const noexcept;
    boost::optional<Path> data_profile() const;
    
private:
    struct Components
    {
        Components() = delete;
        
        Components(ReferenceGenome&& reference, ReadManager&& read_manager,
                   VcfWriter&& output, const options::OptionMap& options);
        
        Components(const Components&)            = delete;
        Components& operator=(const Components&) = delete;
        Components(Components&&)                 = default;
        Components& operator=(Components&&)      = default;
        
        ~Components() = default;
        
        ReferenceGenome reference;
        ReadManager read_manager;
        std::vector<SampleName> samples;
        InputRegionMap regions;
        std::vector<GenomicRegion::ContigName> contigs;
        boost::optional<ReadSetProfile> reads_profile;
        ReadPipe read_pipe;
        HaplotypeLikelihoodModel haplotype_likelihood_model;
        CallerFactory caller_factory;
        boost::optional<ReadPipe> filter_read_pipe;
        VcfWriter output;
        boost::optional<VcfWriter> filtered_output;
        boost::optional<unsigned> num_threads;
        MemoryFootprint read_buffer_footprint;
        std::size_t read_buffer_size;
        ProgressMeter progress_meter;
        PloidyMap ploidies;
        boost::optional<Pedigree> pedigree;
        bool sites_only;
        boost::optional<Path> legacy;
        boost::optional<Path> filter_request;
        boost::optional<Path> bamout;
        BAMRealigner::Config bamout_config;
        boost::optional<Path> data_profile;
        // Components that require temporary directory during construction appear last to make
        // exception handling easier.
        boost::optional<Path> temp_directory;
        std::unique_ptr<VariantCallFilterFactory> call_filter_factory;
        
        void setup_progress_meter(const options::OptionMap& options);
        void set_read_buffer_size(const options::OptionMap& options);
        void setup_writers(const options::OptionMap& options);
        void setup_filter_read_pipe(const options::OptionMap& options);
    };
    
    Components components_;
    
    void update_dependents() noexcept;
};

GenomeCallingComponents collate_genome_calling_components(const options::OptionMap& options);

bool validate(const GenomeCallingComponents& components);

void cleanup(GenomeCallingComponents& components) noexcept;

struct ContigCallingComponents
{
    std::reference_wrapper<const ReferenceGenome> reference;
    std::reference_wrapper<const ReadManager> read_manager;
    InputRegionMap::mapped_type regions;
    std::reference_wrapper<const std::vector<SampleName>> samples;
    std::unique_ptr<const Caller> caller;
    std::size_t read_buffer_size;
    std::reference_wrapper<VcfWriter> output;
    std::reference_wrapper<ProgressMeter> progress_meter;
    
    ContigCallingComponents() = delete;
    
    ContigCallingComponents(const GenomicRegion::ContigName& contig,
                            GenomeCallingComponents& genome_components);
    
    ContigCallingComponents(const GenomicRegion::ContigName& contig, VcfWriter& output,
                            GenomeCallingComponents& genome_components);
    
    ContigCallingComponents(const ContigCallingComponents&)            = delete;
    ContigCallingComponents& operator=(const ContigCallingComponents&) = delete;
    ContigCallingComponents(ContigCallingComponents&&)                 = default;
    ContigCallingComponents& operator=(ContigCallingComponents&&)      = default;
    
    ~ContigCallingComponents() = default;
};

} // namespace octopus

#endif
