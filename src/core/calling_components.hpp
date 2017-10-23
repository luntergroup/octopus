// Copyright (c) 2016 Daniel Cooke
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
#include "io/reference/reference_genome.hpp"
#include "io/read/read_manager.hpp"
#include "io/variant/vcf_writer.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "core/callers/caller_factory.hpp"
#include "core/csr/filters/variant_call_filter_factory.hpp"
#include "core/csr/filters/variant_call_filter.hpp"
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
    std::size_t read_buffer_size() const noexcept;
    const boost::optional<Path>& temp_directory() const noexcept;
    boost::optional<unsigned> num_threads() const noexcept;
    const CallerFactory& caller_factory() const noexcept;
    boost::optional<VcfWriter&> filtered_output();
    std::unique_ptr<csr::VariantCallFilter> call_filter() const noexcept;
    ProgressMeter& progress_meter() noexcept;
    bool sites_only() const noexcept;
    boost::optional<Path> legacy() const;
    
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
        ReadPipe read_pipe;
        CallerFactory caller_factory;
        boost::optional<VariantCallFilterFactory> call_filter_factory;
        VcfWriter output;
        boost::optional<unsigned> num_threads;
        std::size_t read_buffer_size;
        boost::optional<Path> temp_directory;
        ProgressMeter progress_meter;
        bool sites_only;
        boost::optional<VcfWriter> filtered_output;
        boost::optional<Path> legacy;
        
        void setup_progress_meter(const options::OptionMap& options);
        void set_read_buffer_size(const options::OptionMap& options);
        void setup_writers(const options::OptionMap& options);
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
    const InputRegionMap::mapped_type regions;
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
