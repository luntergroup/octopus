// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "calling_components.hpp"

#include <sstream>
#include <utility>
#include <iterator>
#include <algorithm>
#include <functional>
#include <exception>

#include "config/config.hpp"
#include "config/option_collation.hpp"
#include "utils/map_utils.hpp"
#include "utils/thread_pool.hpp"
#include "logging/logging.hpp"
#include "exceptions/user_error.hpp"

namespace octopus {

namespace {

namespace fs = boost::filesystem;

VcfWriter make_vcf_writer(boost::optional<fs::path> dst)
{
    return dst ? VcfWriter {std::move(*dst)} : VcfWriter {};
}

} // namespace

GenomeCallingComponents::GenomeCallingComponents(ReferenceGenome&& reference, ReadManager&& read_manager,
                                                 VcfWriter&& output, const options::OptionMap& options)
: components_ {std::move(reference), std::move(read_manager), std::move(output), options}
{}

GenomeCallingComponents::GenomeCallingComponents(GenomeCallingComponents&& other) noexcept
:  components_ {std::move(other.components_)}
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

MemoryFootprint GenomeCallingComponents::read_buffer_footprint() const noexcept
{
    return components_.read_buffer_footprint;
}

std::size_t GenomeCallingComponents::read_buffer_size() const noexcept
{
    return components_.read_buffer_size;
}

const boost::optional<GenomeCallingComponents::Path>& GenomeCallingComponents::temp_directory() const noexcept
{
    return components_.temp_directory;
}

bool GenomeCallingComponents::keep_temporary_files() const noexcept
{
    return components_.keep_temp_files;
}

boost::optional<unsigned> GenomeCallingComponents::num_threads() const noexcept
{
    return components_.num_threads;
}

const HaplotypeLikelihoodModel& GenomeCallingComponents::haplotype_likelihood_model() const noexcept
{
    return components_.haplotype_likelihood_model;
}

HaplotypeLikelihoodModel GenomeCallingComponents::realignment_haplotype_likelihood_model() const
{
    return components_.realignment_haplotype_likelihood_model;
}

const CallerFactory& GenomeCallingComponents::caller_factory() const noexcept
{
    return components_.caller_factory;
}

boost::optional<VcfWriter&> GenomeCallingComponents::filtered_output() noexcept
{
    if (components_.filtered_output) {
        return *components_.filtered_output; // convert to reference
    } else {
        return boost::none;
    }
}

boost::optional<const VcfWriter&> GenomeCallingComponents::filtered_output() const noexcept
{
    if (components_.filtered_output) {
        return *components_.filtered_output; // convert to reference
    } else {
        return boost::none;
    }
}

const VariantCallFilterFactory& GenomeCallingComponents::call_filter_factory() const
{
    return *components_.call_filter_factory;
}

ReadPipe& GenomeCallingComponents::filter_read_pipe() noexcept
{
    return components_.filter_read_pipe ? *components_.filter_read_pipe : read_pipe();
}

const ReadPipe& GenomeCallingComponents::filter_read_pipe() const noexcept
{
    return components_.filter_read_pipe ? *components_.filter_read_pipe : read_pipe();
}

ProgressMeter& GenomeCallingComponents::progress_meter() noexcept
{
    return components_.progress_meter;
}

boost::optional<GenomeCallingComponents::Path> GenomeCallingComponents::filter_request() const
{
    return components_.filter_request;
}

boost::optional<GenomeCallingComponents::Path> GenomeCallingComponents::bamout() const
{
    return components_.bamout;
}

BAMRealigner::Config GenomeCallingComponents::bamout_config() const noexcept
{
    return components_.bamout_config;
}

boost::optional<const ReadSetProfile&> GenomeCallingComponents::reads_profile() const noexcept
{
    if (components_.reads_profile) {
        return *components_.reads_profile;
    } else {
        return boost::none;
    }
}

boost::optional<GenomeCallingComponents::Path> GenomeCallingComponents::data_profile() const
{
    return components_.data_profile;
}

IndelProfiler::ProfileConfig GenomeCallingComponents::profiler_config() const
{
    return components_.profiler_config;
}

bool GenomeCallingComponents::sites_only() const noexcept
{
    return components_.sites_only;
}

const PloidyMap& GenomeCallingComponents::ploidies() const noexcept
{
    return components_.ploidies;
}

boost::optional<Pedigree> GenomeCallingComponents::pedigree() const
{
    return components_.pedigree;
}

namespace {

std::vector<ContigName>
copy_unmapped_contigs(const ReadManager& rm, const ReferenceGenome& reference, const std::vector<ContigName>& contigs)
{
    std::vector<ContigName> result{};
    result.reserve(contigs.size());
    std::copy_if(std::cbegin(contigs), std::cend(contigs), std::back_inserter(result),
                 [&](const auto& contig) {
                     try {
                         rm.has_reads(reference.contig_region(contig));
                         return false;
                     } catch (...) {
                         return true;
                     }
                 });
    return result;
}

std::vector<ContigName> get_unmapped_contigs(const ReadManager& rm, const ReferenceGenome& reference)
{
    return copy_unmapped_contigs(rm, reference, reference.contig_names());
}

auto get_search_regions(const options::OptionMap& options, const ReferenceGenome& reference, const ReadManager& rm)
{
    auto result = options::get_search_regions(options, reference);
    if (options::ignore_unmapped_contigs(options)) {
        const auto unmapped_contigs = get_unmapped_contigs(rm, reference);
        if (!unmapped_contigs.empty()) {
            logging::WarningLogger warn_log{};
            stream(warn_log) << "Ignoring " << unmapped_contigs.size() << " unmapped contigs";
            for (const auto& contig : unmapped_contigs) {
                result.erase(contig);
            }
        }
    }
    return result;
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
        case ContigOutputOrder::lexicographicalAscending:
            result = std::less<> {};
            break;
        case ContigOutputOrder::lexicographicalDescending:
            result = std::greater<> {};
            break;
        case ContigOutputOrder::contigSizeAscending:
            result = [&reference](const auto& lhs, const auto& rhs) -> bool {
                return reference.contig_size(lhs) < reference.contig_size(rhs);
            };
            break;
        case ContigOutputOrder::contigSizeDescending:
            result = [&reference](const auto& lhs, const auto& rhs) -> bool {
                return reference.contig_size(lhs) > reference.contig_size(rhs);
            };
            break;
        case ContigOutputOrder::referenceIndex: {
            auto reference_contigs = reference.contig_names();
            result = [reference_contigs = std::move(reference_contigs)]
            (const auto& lhs, const auto& rhs) -> bool {
                return index_of(reference_contigs, lhs) < index_of(reference_contigs, rhs);
            };
            break;
        }
        case ContigOutputOrder::referenceIndexReversed: {
            auto reference_contigs = reference.contig_names();
            result = [reference_contigs = std::move(reference_contigs)]
            (const auto& lhs, const auto& rhs) -> bool {
                return index_of(reference_contigs, lhs) < index_of(reference_contigs, rhs);
            };
            break;
        }
        case ContigOutputOrder::unspecified:
            result = std::less<> {};
            break;
    }
    
    return result;
}

auto get_contigs(const InputRegionMap& regions, const ReferenceGenome& reference,
                 const options::ContigOutputOrder order)
{
    auto result = extract_keys(regions);
    const auto cmp = get_sorter(order, reference);
    // Pass a lambda here rather than cmp to avoid compiler bug
    // https://llvm.org/bugs/show_bug.cgi?id=24115
    std::sort(std::begin(result), std::end(result),
              [&cmp](const auto& lhs, const auto& rhs) {
                  return cmp(lhs, rhs);
              });
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
                                       [&file_samples](const auto& sample) {
                                           return is_in_file_samples(sample, file_samples);
                                       });
        const auto num_not_found = std::distance(it, std::end(*user_samples));
        if (num_not_found > 0) {
            std::ostringstream ss{};
            ss << "The requested calling sample";
            if (num_not_found > 1) ss << 's';
            ss << " ";
            std::transform(it, std::end(*user_samples), std::ostream_iterator<SampleName>(ss, ", "),
                           [](auto sample) { return "'" + sample + "'"; });
            if (num_not_found == 1) {
                ss << "is";
            } else {
                ss << "are";
            }
            ss << " not present in any of the read files";
            logging::WarningLogger log{};
            log << ss.str();
            user_samples->erase(it, std::end(*user_samples));
        }
        return *user_samples;
    } else {
        return file_samples;
    }
}

void drop_unused_samples(std::vector<SampleName>& calling_samples, ReadManager& rm)
{
    if (rm.num_samples() <= calling_samples.size()) return;
    std::sort(std::begin(calling_samples), std::end(calling_samples));
    auto managed_samples = rm.samples();
    assert(calling_samples.size() < managed_samples.size());
    std::sort(std::begin(managed_samples), std::end(managed_samples));
    std::vector<SampleName> unused_samples{};
    unused_samples.reserve(managed_samples.size() - calling_samples.size());
    std::set_difference(std::cbegin(managed_samples), std::cend(managed_samples),
                        std::cbegin(calling_samples), std::cend(calling_samples),
                        std::back_inserter(unused_samples));
    rm.drop_samples(unused_samples);
}

bool is_multithreaded_run(const options::OptionMap& options) noexcept
{
    const auto num_threads = options::get_num_threads(options);
    return !num_threads || *num_threads > 1;
}

bool is_stdout_output(const options::OptionMap& options)
{
    return !options::get_output_path(options);
}

bool require_temp_dir_for_filtering(const options::OptionMap& options)
{
    return options::is_call_filtering_requested(options)
           && (!options::keep_unfiltered_calls(options) || is_stdout_output(options));
}

bool is_temp_directory_needed(const options::OptionMap& options)
{
    return is_multithreaded_run(options) || require_temp_dir_for_filtering(options);
}

boost::optional<fs::path> get_temp_directory(const options::OptionMap& options)
{
    if (is_temp_directory_needed(options)) {
        return options::create_temp_file_directory(options);
    } else {
        return boost::none;
    }
}

namespace {

bool includes(const std::vector<unsigned>& values, const unsigned value) noexcept
{
    return std::find(std::cbegin(values), std::cend(values), value) != std::cend(values);
}

auto profile_reads_helper(const std::vector<SampleName>& samples,
                          const ReferenceGenome& reference,
                          const InputRegionMap& input_regions,
                          const ReadManager& source,
                          const PloidyMap& ploidies,
                          const options::OptionMap& options)
{
    ReadSetProfileConfig config {};
    config.fragment_size = options::max_read_length(options);
    const auto threads = options::get_num_threads(options);
    ThreadPool pool {threads ? *threads : 1};
    boost::optional<ThreadPool&> workers {};
    if (threads && *threads > 1) {
        workers = pool;
    }
    if (samples.size() == 1) {
        auto result = profile_reads(samples, reference, input_regions, source, config, workers);
        if (result) result->depth_stats.sample.clear(); // no need to keep this duplicate info
        return result;
    } else if (options::use_same_read_profile_for_all_samples(options)) {
        std::vector<SampleName> profile_samples {};
        std::unordered_map<GenomicRegion::ContigName, std::vector<unsigned>> ploidies_covered {};
        for (const auto& sample : samples) {
            bool include_sample {false};
            for (const auto& contig : extract_keys(input_regions)) {
                const auto ploidy = ploidies.of(sample, contig);
                if (ploidy > 0 && !includes(ploidies_covered[contig], ploidy)) {
                    ploidies_covered[contig].push_back(ploidy);
                    include_sample = true;
                }
            }
            if (include_sample) profile_samples.push_back(sample);
        }
        auto result = profile_reads(profile_samples, reference, input_regions, source, config, workers);
        if (result) result->depth_stats.sample.clear();
        return result;
    } else {
        return profile_reads(samples, reference, input_regions, source, config, workers);
    }
}

static const AlignedRead typical_illumina_read {
    "HISEQ1:9:H8962ADXX:2:1108:11915:94551",
    GenomicRegion {"1", 63492953, 63493103},
    "GGCCAGAGAGAGAGTAGGTGAATCTGATCTCAGAATGTAAGCTCCTGACCAGTACAGCAGCCTGCATGCCCCCAGGAGCTGGCAGAGGAGGAGGAGGAGGAGGAGGCGGCAGATGATTCACAGCCACAATAGCATTGGCAACACTGGG",
    AlignedRead::BaseQualityVector {33,35,35,35,34,32,27,35,35,35,35,35,36,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,36,35,36,36,36,35,35,35,35,35,34,34,32,35,35,34,34,34,34,35,35,35,35,35,34,34,36,36,36,34,36,34,35,37,37,37,39,39,39,38,39,39,41,41,41,41,41,41,41,41,41,41,41,41,40,41,41,41,41,41,41,41,41,40,38,40,40,40,41,40,40,40,40,41,41,41,40,41,41,41,41,40,40,40,41,41,41,41,41,41,41,41,41,41,41,41,40,39,41,41,41,41,41,41,41,41,39,39,39,39,39,37,37,37,37,37,34,34,34},
    parse_cigar("83M3D65M"),
    60,
    AlignedRead::Flags {},
    "None",
    std::vector<std::pair<AlignedRead::Tag, AlignedRead::Annotation>> {}
};

} // namespace

auto estimate_read_memory_footprint(const boost::optional<ReadSetProfile>& profile) noexcept
{
    MemoryFootprint result;
    if (profile) {
        result = std::min(profile->memory_stats.mean + profile->memory_stats.stdev, profile->memory_stats.max);
        if (profile->fragmented_memory_stats) {
            result += profile->fragmented_memory_stats->median;
        }
    } else {
        result = footprint(typical_illumina_read);
        logging::WarningLogger log {};
        log << "Could not estimate read size from data, resorting to default";
    }
    auto debug_log = logging::get_debug_log();
    if (debug_log) stream(*debug_log) << "Estimated read memory footprint is " << result;
    return result;
}

std::size_t calculate_max_num_reads(MemoryFootprint max_buffer_size, const boost::optional<ReadSetProfile>& profile) noexcept
{
    static constexpr MemoryFootprint min_buffer_size {50'000'000}; // 50Mb
    if (max_buffer_size < min_buffer_size) {
        static bool warned {false};
        if (!warned) {
            logging::WarningLogger warn_log {};
            stream(warn_log) << "Ignoring given maximum read buffer size of " << max_buffer_size
                             << " as this size is too small. Setting maximum to "
                             << min_buffer_size << " instead.";
            warned = true;
        }
        max_buffer_size = min_buffer_size;
    }
    auto estimated_read_footprint = estimate_read_memory_footprint(profile);
    return max_buffer_size.bytes() / estimated_read_footprint.bytes();
}

auto add_identifier(const fs::path& base, const std::string& identifier)
{
    const auto old_stem  = base.stem();
    const auto extension = base.extension();
    fs::path new_stem;
    if (extension.string() == ".gz") {
        new_stem = old_stem.stem().string() + "." + identifier + old_stem.extension().string() + extension.string();
    } else {
        new_stem = old_stem.string() + "." + identifier + extension.string();
    }
    return base.parent_path() / new_stem;
}

auto get_unfiltered_path(const fs::path& native)
{
    return add_identifier(native, "unfiltered");
}

auto generate_temp_output_path(const fs::path& temp_directory)
{
    return temp_directory / "octopus_unfiltered.bcf";
}

bool all_samples_in_vcf(std::vector<SampleName> samples, const VcfReader& in)
{
    std::sort(std::begin(samples), std::end(samples));
    auto vcf_samples = in.fetch_header().samples();
    std::sort(std::begin(vcf_samples), std::end(vcf_samples));
    return samples == vcf_samples;
}

bool all_samples_in_vcf(const std::vector<SampleName>& samples, const fs::path& vcf_in)
{
    const VcfReader tmp {vcf_in};
    return all_samples_in_vcf(samples, tmp);
}

class InputVCFError : public UserError
{
    std::string do_where() const override { return "GenomeCallingComponents"; }
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "the VCF file you specified ";
        ss << vcf_path_;
        ss << " contains samples not in the input BAM files";
        return ss.str();
    }
    std::string do_help() const override
    {
        return "ensure the input VCF file samples match the input BAM samples";
    }
    
    fs::path vcf_path_;

public:
    InputVCFError(fs::path vcf_path) : vcf_path_ {std::move(vcf_path)} {}
    
    virtual ~InputVCFError() override = default;
};

template <typename T>
boost::optional<const T&> optional_cref(const boost::optional<T>& v)
{
    if (v) {
        return *v;
    } else {
        return boost::none;
    }
}

} // namespace

GenomeCallingComponents::Components::Components(ReferenceGenome&& reference, ReadManager&& read_manager,
                                                VcfWriter&& output, const options::OptionMap& options)
: reference {std::move(reference)}
, read_manager {std::move(read_manager)}
, samples {extract_samples(options, this->read_manager)}
, regions {get_search_regions(options, this->reference, this->read_manager)}
, contigs {get_contigs(this->regions, this->reference, options::get_contig_output_order(options))}
, ploidies {options::get_ploidy_map(options)}
, reads_profile {profile_reads_helper(this->samples, this->reference, this->regions, this->read_manager, this->ploidies, options)}
, read_pipe {options::make_read_pipe(this->read_manager, this->reference, this->samples, options)}
, haplotype_likelihood_model {options::make_calling_haplotype_likelihood_model(options, optional_cref(this->reads_profile))}
, realignment_haplotype_likelihood_model {options::make_realignment_haplotype_likelihood_model(haplotype_likelihood_model, optional_cref(this->reads_profile), options)}
, caller_factory {options::make_caller_factory(this->reference, this->read_pipe, this->regions, options, optional_cref(this->reads_profile))}
, filter_read_pipe {}
, output {std::move(output)}
, filtered_output {}
, num_threads {options::get_num_threads(options)}
, read_buffer_footprint {options::get_target_read_buffer_size(options)}
, read_buffer_size {}
, progress_meter {regions}
, pedigree {options::get_pedigree(options, samples)}
, sites_only {options::call_sites_only(options)}
, filter_request {}
, bamout {options::bamout_request(options)}
, bamout_config {}
, data_profile {options::data_profile_request(options)}
, profiler_config {}
{
    drop_unused_samples(this->samples, this->read_manager);
    setup_progress_meter(options);
    set_read_buffer_size(options);
    setup_filter_read_pipe(options);
    filter_request = options::filter_request(options);
    if (filter_request && !all_samples_in_vcf(samples, *filter_request)) {
        throw InputVCFError {*filter_request};
    }
    temp_directory = get_temp_directory(options);
    try {
        call_filter_factory = options::make_call_filter_factory(this->reference, this->read_pipe, options, this->temp_directory);
        setup_writers(options);
    } catch (...) {
        if (temp_directory) fs::remove_all(*temp_directory);
        throw;
    }
    keep_temp_files = options::keep_temporary_files(options);
    bamout_config.alignment_model = realignment_haplotype_likelihood_model;
    bamout_config.copy_hom_ref_reads = options::full_bamouts_requested(options);
    bamout_config.max_buffer = read_buffer_footprint;
    bamout_config.max_threads = num_threads;
    bamout_config.read_linkage.linkage = options::get_read_linkage_type(options);
    profiler_config.alignment_model = bamout_config.alignment_model;
    if (reads_profile && reads_profile->length_stats.median > 1'000) {
        profiler_config.ignore_likely_misaligned_reads = false;
    }
}

void GenomeCallingComponents::Components::setup_progress_meter(const options::OptionMap& options)
{
    const auto num_bp_to_process = sum_region_sizes(regions);
    if (num_bp_to_process < 100000000) {
        progress_meter.set_max_tick_size(1.0);
    } else if (num_bp_to_process < 1000000000) {
        progress_meter.set_max_tick_size(0.5);
    } else {
        progress_meter.set_max_tick_size(0.1);
    }
}

void GenomeCallingComponents::Components::set_read_buffer_size(const options::OptionMap& options)
{
    if (!samples.empty() && !regions.empty() && read_manager.good()) {
        read_buffer_size = calculate_max_num_reads(options::get_target_read_buffer_size(options), reads_profile);
    }
}

void GenomeCallingComponents::Components::setup_writers(const options::OptionMap& options)
{
    if (call_filter_factory) {
        const auto final_output_path = output.path();
        filtered_output = std::move(output);
        fs::path prefilter_path;
        if (final_output_path) {
            if (options::keep_unfiltered_calls(options)) {
                prefilter_path = get_unfiltered_path(*final_output_path);
            } else {
                assert(temp_directory);
                prefilter_path = *temp_directory / get_unfiltered_path(final_output_path->filename());
            }
        } else {
            assert(temp_directory);
            prefilter_path = generate_temp_output_path(*temp_directory);
        }
        output.open(std::move(prefilter_path));
    }
}

void GenomeCallingComponents::Components::setup_filter_read_pipe(const options::OptionMap& options)
{
    if (!options::use_calling_read_pipe_for_call_filtering(options)) {
        filter_read_pipe = options::make_call_filter_read_pipe(read_manager, reference, samples, options);
    }
}

void GenomeCallingComponents::update_dependents() noexcept
{
    components_.read_pipe.set_read_manager(components_.read_manager);
    if (components_.filter_read_pipe) {
        components_.filter_read_pipe->set_read_manager(components_.read_manager);
    }
    components_.caller_factory.set_reference(components_.reference);
    components_.caller_factory.set_read_pipe(components_.read_pipe);
}

namespace {

bool all_reference_contigs_mapped(const ReadManager& rm, const ReferenceGenome& reference)
{
    const auto contigs = reference.contig_names();
    return std::all_of(std::cbegin(contigs), std::cend(contigs),
                       [&reference, &rm] (const auto& contig) {
                           bool ok {true};
                           try {
                               rm.has_reads(reference.contig_region(contig));
                           } catch (...) {
                               ok = false;
                           }
                           return ok;
                       });
}

class UnmatchedReference : public UserError
{
public:
    UnmatchedReference(const ReferenceGenome& reference)
    : reference_name_ {reference.name()}
    , why_ {"Some or all of the contigs in the reference genome (" +  reference_name_
            + ") are not present in the read files"}
    {}

private:
    std::string do_where() const override
    {
        return "validate";
    }
    
    std::string do_why() const override
    {
        return why_;
    }
    
    std::string do_help() const override
    {
        return "Ensure the reference genome used for mapping is the same as the one used for calling"
        " (" + reference_name_ + ") and all input contigs are present in the read headers";
    }
    
    std::string reference_name_, why_;
};

VcfWriter make_output_vcf_writer(const options::OptionMap& options)
{
    return make_vcf_writer(options::get_output_path(options));
}

} // namespace

GenomeCallingComponents collate_genome_calling_components(const options::OptionMap& options)
{
    auto reference    = options::make_reference(options);
    auto read_manager = options::make_read_manager(options);
    // Check this here to avoid creating output file on error
    if (!options::ignore_unmapped_contigs(options) && !all_reference_contigs_mapped(read_manager, reference)) {
        throw UnmatchedReference {reference};
    }
    auto output = make_output_vcf_writer(options);
    if (output.path() && output.path()->extension().string() == ".vcf") {
        logging::WarningLogger log {};
        log << "Uncompressed VCF output requested - this may result in large output files";
    }
    return GenomeCallingComponents {
        std::move(reference),
        std::move(read_manager),
        std::move(output),
        options
    };
}

bool validate(const GenomeCallingComponents& components)
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
    if (components.temp_directory() && !components.keep_temporary_files()) {
        try {
            const auto num_files_removed = fs::remove_all(*components.temp_directory());
            stream(log) << "Removed " << num_files_removed << " temporary files";
        } catch (const std::exception& e) {
            stream(log) << "Cleanup failed with exception: " << e.what();
        }
    }
}

ContigCallingComponents::ContigCallingComponents(const GenomicRegion::ContigName& contig,
                                                 GenomeCallingComponents& genome_components)
: reference {genome_components.reference()}
, read_manager {genome_components.read_manager()}
, regions {genome_components.search_regions().at(contig)}
, samples {genome_components.samples()}
, caller {genome_components.caller_factory().make(contig)}
, read_buffer_size {genome_components.read_buffer_size()}
, output {genome_components.output()}
, progress_meter {genome_components.progress_meter()}
{}

ContigCallingComponents::ContigCallingComponents(const GenomicRegion::ContigName& contig, VcfWriter& output,
                                                 GenomeCallingComponents& genome_components)
: reference {genome_components.reference()}
, read_manager {genome_components.read_manager()}
, regions {genome_components.search_regions().at(contig)}
, samples {genome_components.samples()}
, caller {genome_components.caller_factory().make(contig)}
, read_buffer_size {genome_components.read_buffer_size()}
, output {output}
, progress_meter {genome_components.progress_meter()}
{}

} // namespace octopus
