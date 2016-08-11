// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "calling_components.hpp"

#include <sstream>
#include <utility>
#include <iterator>
#include <algorithm>
#include <functional>
#include <exception>

#include <config/option_collation.hpp>
#include <utils/read_size_estimator.hpp>
#include <utils/map_utils.hpp>
#include <logging/logging.hpp>

namespace octopus {

GenomeCallingComponents::GenomeCallingComponents(ReferenceGenome&& reference, ReadManager&& read_manager,
                        VcfWriter&& output, const options::OptionMap& options)
:
components_ {std::move(reference), std::move(read_manager), std::move(output), options}
{}

GenomeCallingComponents::GenomeCallingComponents(GenomeCallingComponents&& other) noexcept
:
components_ {std::move(other.components_)}
{
    update_dependents();
}

const ReferenceGenome& GenomeCallingComponents::reference() const noexcept
{
    return components_.reference;
}

ReadManager& GenomeCallingComponents::read_manager() noexcept
{
    return components_.read_manager;
}

const ReadManager& GenomeCallingComponents::read_manager() const noexcept
{
    return components_.read_manager;
}

ReadPipe& GenomeCallingComponents::read_pipe() noexcept
{
    return components_.read_pipe;
}

const ReadPipe& GenomeCallingComponents::read_pipe() const noexcept
{
    return components_.read_pipe;
}

const std::vector<SampleName>& GenomeCallingComponents::samples() const noexcept
{
    return components_.samples;
}

const InputRegionMap& GenomeCallingComponents::search_regions() const noexcept
{
    return components_.regions;
}

const std::vector<GenomicRegion::ContigName>& GenomeCallingComponents::contigs() const noexcept
{
    return components_.contigs;
}

VcfWriter& GenomeCallingComponents::output() noexcept
{
    return components_.output;
}

const VcfWriter& GenomeCallingComponents::output() const noexcept
{
    return components_.output;
}

std::size_t GenomeCallingComponents::read_buffer_size() const noexcept
{
    return components_.read_buffer_size;
}

const boost::optional<boost::filesystem::path>& GenomeCallingComponents::temp_directory() const noexcept
{
    return components_.temp_directory;
}

boost::optional<unsigned> GenomeCallingComponents::num_threads() const noexcept
{
    return components_.num_threads;
}

const CallerFactory& GenomeCallingComponents::caller_factory() const noexcept
{
    return components_.caller_factory;
}

ProgressMeter& GenomeCallingComponents::progress_meter() noexcept
{
    return components_.progress_meter;
}

template <typename T>
std::size_t index_of(const std::vector<T>& elements, const T& value)
{
    return std::distance(std::cbegin(elements), std::find(std::cbegin(elements), std::cend(elements), value));
}
    
auto get_sorter(const options::ContigOutputOrder order, const ReferenceGenome& reference)
{
    using options::ContigOutputOrder;
    
    std::function<bool(const ContigName&, const ContigName&)> result;
    
    switch (order) {
        case ContigOutputOrder::LexicographicalAscending:
            result = std::less<> {};
            break;
        case ContigOutputOrder::LexicographicalDescending:
            result = std::greater<> {};
            break;
        case ContigOutputOrder::ContigSizeAscending:
            result = [&reference] (const auto& lhs, const auto& rhs) -> bool {
                return reference.contig_size(lhs) < reference.contig_size(rhs);
            };
            break;
        case ContigOutputOrder::ContigSizeDescending:
            result = [&reference] (const auto& lhs, const auto& rhs) -> bool {
                return reference.contig_size(lhs) > reference.contig_size(rhs);
            };
            break;
        case ContigOutputOrder::AsInReferenceIndex:
        {
            auto reference_contigs = reference.contig_names();
            result = [&reference, reference_contigs = std::move(reference_contigs)]
            (const auto& lhs, const auto& rhs) -> bool {
                return index_of(reference_contigs, lhs) < index_of(reference_contigs, rhs);
            };
            break;
        }
        case ContigOutputOrder::AsInReferenceIndexReversed:
        {
            auto reference_contigs = reference.contig_names();
            result = [&reference, reference_contigs = std::move(reference_contigs)]
            (const auto& lhs, const auto& rhs) -> bool {
                return index_of(reference_contigs, lhs) < index_of(reference_contigs, rhs);
            };
            break;
        }
        case ContigOutputOrder::Unspecified:
            result = std::less<> {};
            break;
    }
    
    return result;
}

auto get_contigs(const InputRegionMap& regions, const ReferenceGenome& reference,
                 const options::ContigOutputOrder order)
{
    auto result = extract_keys(regions);
    
    std::sort(std::begin(result), std::end(result), get_sorter(order, reference));
    
    return result;
}

template <typename Container>
bool is_in_file_samples(const SampleName& sample, const Container& file_samples)
{
    return std::find(std::cbegin(file_samples), std::cend(file_samples), sample) != std::cend(file_samples);
}

std::vector<SampleName> extract_samples(const options::OptionMap& options, const ReadManager& read_manager)
{
    auto user_samples = options::get_user_samples(options);
    auto file_samples = read_manager.samples();
    
    if (user_samples) {
        const auto it = std::partition(std::begin(*user_samples), std::end(*user_samples),
                                       [&file_samples] (const auto& sample) {
                                           return is_in_file_samples(sample, file_samples);
                                       });
        
        const auto num_not_found = std::distance(it, std::end(*user_samples));
        
        if (num_not_found > 0) {
            std::ostringstream ss {};
            ss << "The requested calling sample";
            if (num_not_found > 1) ss << 's';
            ss << " ";
            std::transform(it, std::end(*user_samples), std::ostream_iterator<SampleName>(ss, ", "),
                           [] (auto sample) { return "'" + sample + "'"; });
            if (num_not_found == 1) {
                ss << "is";
            } else {
                ss << "are";
            }
            ss << " not present in any of the read files";
            logging::WarningLogger log {};
            log << ss.str();
            user_samples->erase(it, std::end(*user_samples));
        }
        
        return *user_samples;
    } else {
        return file_samples;
    }
}

auto estimate_read_size(const std::vector<SampleName>& samples,
                        const InputRegionMap& input_regions,
                        ReadManager& read_manager)
{
    auto result = estimate_mean_read_size(samples, input_regions, read_manager);
    
    if (!result) {
        logging::WarningLogger log {};
        log << "Could not estimate read size from data, resorting to default";
        
        return default_read_size_estimate();
    }
    
    return *result;
}

std::size_t calculate_max_num_reads(const std::size_t max_buffer_bytes,
                                    const std::vector<SampleName>& samples,
                                    const InputRegionMap& input_regions,
                                    ReadManager& read_manager)
{
    if (samples.empty()) return 0;
    
    static constexpr std::size_t min_buffer_bytes {1'000'000};
    
    return std::max(max_buffer_bytes, min_buffer_bytes) / estimate_read_size(samples, input_regions, read_manager);
}

GenomeCallingComponents::Components::Components(ReferenceGenome&& reference, ReadManager&& read_manager,
                                                VcfWriter&& output, const options::OptionMap& options)
:
reference {std::move(reference)},
read_manager {std::move(read_manager)},
samples {extract_samples(options, this->read_manager)},
regions {options::get_search_regions(options, this->reference)},
contigs {get_contigs(this->regions, this->reference, options::get_contig_output_order(options))},
read_pipe {options::make_read_pipe(this->read_manager, this->samples, options)},
caller_factory {options::make_caller_factory(this->reference, this->read_pipe, this->regions, options)},
output {std::move(output)},
num_threads {options::get_num_threads(options)},
read_buffer_size {},
temp_directory {((!num_threads || *num_threads > 1) ? options::create_temp_file_directory(options) : boost::none)},
progress_meter {regions}
{
    const auto num_bp_to_process = sum_region_sizes(regions);
    
    if (num_bp_to_process < 100000000) {
        progress_meter.set_percent_block_size(1.0);
    } else if (num_bp_to_process < 1000000000) {
        progress_meter.set_percent_block_size(0.5);
    } else {
        progress_meter.set_percent_block_size(0.1);
    }
    
    if (!samples.empty() && !regions.empty() && read_manager.good()) {
        read_buffer_size = calculate_max_num_reads(options::get_target_read_buffer_size(options),
                                                   this->samples, this->regions,
                                                   this->read_manager);
    }
}

void GenomeCallingComponents::update_dependents() noexcept
{
    components_.read_pipe.set_read_manager(components_.read_manager);
    components_.caller_factory.set_reference(components_.reference);
    components_.caller_factory.set_read_pipe(components_.read_pipe);
}

bool check_components_valid(const GenomeCallingComponents& components)
{
    if (components.samples().empty()) {
        logging::WarningLogger log {};
        log << "No samples detected - at least one is required for calling";
        return false;
    }
    
    if (components.search_regions().empty()) {
        logging::WarningLogger log {};
        log << "There are no input regions - at least one is required for calling";
        return false;
    }
    
    return true;
}

void cleanup(GenomeCallingComponents& components) noexcept
{
    logging::InfoLogger log {};
    if (components.temp_directory()) {
        try {
            const auto num_files_removed = fs::remove_all(*components.temp_directory());
            stream(log) << "Removed " << num_files_removed << " temporary files";
        } catch (const std::exception& e) {
            stream(log) << "Cleanup failed with exception: " << e.what();
        }
    }
}

boost::optional<GenomeCallingComponents> collate_genome_calling_components(const options::OptionMap& options)
{
    auto reference = options::make_reference(options);
    
    auto read_manager = options::make_read_manager(options);
    
    auto output = options::make_output_vcf_writer(options);
    
    if (!output.is_open()) {
        return boost::none;
    }
    
    GenomeCallingComponents result {
        std::move(reference),
        std::move(read_manager),
        std::move(output),
        options
    };
    
    check_components_valid(result);
    
    return boost::optional<GenomeCallingComponents> {std::move(result)};
}

ContigCallingComponents::ContigCallingComponents(const GenomicRegion::ContigName& contig,
                                                 GenomeCallingComponents& genome_components)
:
reference {genome_components.reference()},
read_manager {genome_components.read_manager()},
regions {genome_components.search_regions().at(contig)},
samples {genome_components.samples()},
caller {genome_components.caller_factory().make(contig)},
read_buffer_size {genome_components.read_buffer_size()},
output {genome_components.output()},
progress_meter {genome_components.progress_meter()}
{}

ContigCallingComponents::ContigCallingComponents(const GenomicRegion::ContigName& contig, VcfWriter& output,
                                                 GenomeCallingComponents& genome_components)
:
reference {genome_components.reference()},
read_manager {genome_components.read_manager()},
regions {genome_components.search_regions().at(contig)},
samples {genome_components.samples()},
caller {genome_components.caller_factory().make(contig)},
read_buffer_size {genome_components.read_buffer_size()},
output {output},
progress_meter {genome_components.progress_meter()}
{}

} // namespace octopus
