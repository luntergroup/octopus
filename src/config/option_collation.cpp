// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "option_collation.hpp"

#include <string>
#include <iostream>
#include <cctype>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <functional>
#include <utility>
#include <thread>
#include <sstream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>

#include <basics/phred.hpp>
#include <basics/genomic_region.hpp>
#include <basics/aligned_read.hpp>
#include <readpipe/read_pipe_fwd.hpp>
#include <core/tools/coretools.hpp>
#include <core/callers/caller_builder.hpp>
#include <utils/read_stats.hpp>
#include <utils/mappable_algorithms.hpp>
#include <utils/string_utils.hpp>
#include <utils/append.hpp>
#include <utils/maths.hpp>
#include <logging/logging.hpp>
#include <io/region/region_parser.hpp>
#include <io/variant/vcf_reader.hpp>
#include <io/variant/vcf_writer.hpp>
#include <exceptions/user_error.hpp>
#include <exceptions/system_error.hpp>
#include <exceptions/missing_file_error.hpp>

namespace octopus { namespace options {

// unsigned are banned from the option map to prevent user input errors, but once the option
// map is passed they are all safe
unsigned as_unsigned(const std::string& option, const OptionMap& options)
{
    return static_cast<unsigned>(options.at(option).as<int>());
}

bool is_run_command(const OptionMap& options)
{
    return options.count("help") == 0 && options.count("version") == 0;
}

bool is_debug_mode(const OptionMap& options)
{
    return options.at("debug").as<bool>();
}

bool is_trace_mode(const OptionMap& options)
{
    return options.at("trace").as<bool>();
}

boost::optional<fs::path> get_home_directory()
{
    static const auto env_p = std::getenv("HOME");
    
    if (env_p == nullptr) {
        return boost::none;
    }
    
    const fs::path home {env_p};
    
    if (fs::is_directory(home)) {
        return home;
    }
    
    return boost::none;
}

bool is_shorthand_user_path(const fs::path& path)
{
    return !path.empty() && path.string().front() == '~';
}

class UnknownHomeDirectory : public SystemError
{
    std::string do_where() const override
    {
        return "expand_user_path";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "Unable to expand shorthand path you specified ";
        ss << path_;
        ss << " as your home directory cannot be located";
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "ensure your HOME environment variable is set properly";
    }
    
    fs::path path_;
public:
    UnknownHomeDirectory(fs::path p) : path_ {std::move(p)} {}
};

fs::path expand_user_path(const fs::path& path)
{
    if (is_shorthand_user_path(path)) {
        if (path.string().size() > 1 && path.string()[1] == '/') {
            const auto home_dir = get_home_directory();
            if (home_dir) {
                return fs::path {home_dir->string() + path.string().substr(1)};
            }
            throw UnknownHomeDirectory {path};
        }
        return path;
    }
    return path;
}

class InvalidWorkingDirectory : public UserError
{
    std::string do_where() const override
    {
        return "get_working_directory";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The working directory you specified ";
        ss << path_;
        ss << " does not exist";
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "enter a valid working directory";
    }
    
    fs::path path_;
public:
    InvalidWorkingDirectory(fs::path p) : path_ {std::move(p)} {}
};

fs::path get_working_directory(const OptionMap& options)
{
    if (options.count("working-directory") == 1) {
        auto result = expand_user_path(options.at("working-directory").as<std::string>());
        
        if (!fs::exists(result) && !fs::is_directory(result)) {
            throw InvalidWorkingDirectory {result};
        }
        
        return result;
    }
    return fs::current_path();
}

fs::path resolve_path(const fs::path& path, const OptionMap& options)
{
    if (is_shorthand_user_path(path)) {
        return expand_user_path(path); // must be a root path
    }
    
    if (fs::exists(path)) {
        return path; // must be a root path
    }
    
    const auto parent_dir = path.parent_path();
    
    const auto wd = get_working_directory(options);
    
    if (fs::exists(parent_dir) && fs::is_directory(parent_dir)) {
        auto tmp = wd;
        tmp /= path;
        auto wd_parent = tmp.parent_path();
        
        if (fs::exists(wd_parent) && fs::is_directory(wd_parent)) {
            return tmp; // prefer working directory in case of name clash
        }
        
        return path; // must be yet-to-be-created root path
    }
    
    auto result = wd;
    
    result /= path;
    
    return result;
}

namespace {
    struct Line
    {
        std::string line_data;
        
        operator std::string() const
        {
            return line_data;
        }
    };
    
    std::istream& operator>>(std::istream& is, Line& data)
    {
        std::getline(is, data.line_data);
        
        if (!data.line_data.empty() && data.line_data.back() == '\r') {
            data.line_data.pop_back();
        }
        
        return is;
    }
} // namespace

std::vector<fs::path> extract_paths_from_file(const fs::path& file_path, const OptionMap& options)
{
    const auto resolved_path = resolve_path(file_path, options);
    
    std::vector<fs::path> result {};
    
    std::ifstream file {file_path.string()};
    
    if (!file.good()) {
        std::ostringstream ss {};
        ss << "Could not open path file " << file_path;
        throw std::runtime_error {ss.str()};
    }
    
    std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                   std::back_inserter(result), [] (const Line& line) { return line.line_data; });
    
    result.erase(std::remove_if(std::begin(result), std::end(result),
                                [] (const auto& path) { return path.empty(); }),
                 std::end(result));
    
    return result;
}

auto resolve_paths(const std::vector<fs::path>& paths, const OptionMap& options)
{
    std::vector<fs::path> result {};
    result.reserve(paths.size());
    std::transform(std::cbegin(paths), std::cend(paths), std::back_inserter(result),
                   [&] (const auto& path) { return resolve_path(path, options); });
    return result;
}

auto resolve_paths(const std::vector<std::string>& path_strings, const OptionMap& options)
{
    std::vector<fs::path> paths {std::cbegin(path_strings), std::cend(path_strings)};
    return resolve_paths(paths, options);
}

bool is_file_readable(const fs::path& path)
{
    std::ifstream tmp {path.string()};
    return tmp.good();
}

bool is_file_writable(const fs::path& path)
{
    if (!fs::exists(path.parent_path())) {
        return false;
    }
    
    std::ofstream test {path.string()};
    
    const auto result = test.is_open();
    
    fs::remove(path);
    
    return result;
}

bool is_threading_allowed(const OptionMap& options)
{
    unsigned num_threads {1};
    
    if (options.count("threads") == 1) {
        num_threads = as_unsigned("threads", options);
    }
    
    return num_threads != 1;
}

boost::optional<unsigned> get_num_threads(const OptionMap& options)
{
    unsigned num_threads {1};
    
    if (options.count("threads") == 1) {
        num_threads = as_unsigned("threads", options);
    }
    
    if (num_threads > 0) return num_threads;
    
    return boost::none;
}

std::size_t get_target_read_buffer_size(const OptionMap& options)
{
    static constexpr std::size_t scale {1000000000};
    return static_cast<std::size_t>(scale * options.at("target-read-buffer-footprint").as<float>());
}

boost::optional<fs::path> get_debug_log_file_name(const OptionMap& options)
{
    if (options.at("debug").as<bool>()) {
        return resolve_path("octopus_debug.log", options);
    }
    return boost::none;
}

boost::optional<fs::path> get_trace_log_file_name(const OptionMap& options)
{
    if (options.at("trace").as<bool>()) {
        return resolve_path("octopus_trace.log", options);
    }
    return boost::none;
}

ReferenceGenome make_reference(const OptionMap& options)
{
    const fs::path input_path {options.at("reference").as<std::string>()};
    
    auto resolved_path = resolve_path(input_path, options);
    
    const auto ref_cache_size = options.at("max-reference-cache-footprint").as<float>();
    
    static constexpr unsigned Scale {1'000'000};
    
    try {
        return octopus::make_reference(std::move(resolved_path),
                                       static_cast<std::size_t>(Scale * ref_cache_size),
                                       is_threading_allowed(options));
    } catch (MissingFileError& e) {
        e.set_location_specified("the command line option --reference");
        throw;
    } catch (...) {
        throw;
    }
}

InputRegionMap make_search_regions(const std::vector<GenomicRegion>& regions)
{
    InputRegionMap contig_mapped_regions {};
    
    for (const auto& region : regions) {
        contig_mapped_regions[region.contig_name()].insert(region);
    }
    
    InputRegionMap result {};
    result.reserve(contig_mapped_regions.size());
    
    for (const auto& p : contig_mapped_regions) {
        auto covered_contig_regions = extract_covered_regions(p.second);
        result.emplace(std::piecewise_construct,
                       std::forward_as_tuple(p.first),
                       std::forward_as_tuple(std::make_move_iterator(std::begin(covered_contig_regions)),
                                             std::make_move_iterator(std::end(covered_contig_regions))));
    }
    
    return result;
}

InputRegionMap extract_search_regions(const ReferenceGenome& reference)
{
    return make_search_regions(get_all_contig_regions(reference));
}

auto cut(const MappableFlatSet<GenomicRegion>& mappables, const MappableFlatSet<GenomicRegion>& regions)
{
    if (mappables.empty()) return regions;
    
    MappableFlatSet<GenomicRegion> result {};
    
    for (const auto& region : regions) {
        auto overlapped = mappables.overlap_range(region);
        
        if (empty(overlapped)) {
            result.emplace(region);
        } else if (!is_same_region(region, overlapped.front())) {
            auto spliced = region;
            
            if (begins_before(overlapped.front(), spliced)) {
                spliced = right_overhang_region(spliced, overlapped.front());
                overlapped.advance_begin(1);
            }
            
            std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&] (const auto& region) {
                result.emplace(left_overhang_region(spliced, region));
                spliced = expand_lhs(spliced, -begin_distance(spliced, region));
            });
            
            if (ends_before(overlapped.back(), spliced)) {
                result.emplace(right_overhang_region(spliced, overlapped.back()));
            }
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

InputRegionMap extract_search_regions(const std::vector<GenomicRegion>& regions,
                                      std::vector<GenomicRegion>& skip_regions)
{
    auto input_regions = make_search_regions(regions);
    
    const auto skipped = make_search_regions(skip_regions);
    
    InputRegionMap result {input_regions.size()};
    
    for (auto& p : input_regions) {
        if (skipped.count(p.first) == 1) {
            result.emplace(p.first, cut(skipped.at(p.first), std::move(p.second)));
        } else {
            result.emplace(p.first, std::move(p.second));
        }
    }
    
    for (auto it = std::begin(result); it != std::end(result); ) {
        if (it->second.empty()) {
            it = result.erase(it);
        } else {
            ++it;
        }
    }
    
    for (auto& p : result) {
        p.second.shrink_to_fit();
    }
    
    return result;
}

InputRegionMap extract_search_regions(const ReferenceGenome& reference,
                                      std::vector<GenomicRegion>& skip_regions)
{
    return extract_search_regions(get_all_contig_regions(reference), skip_regions);
}
        
std::vector<GenomicRegion> parse_regions(const std::vector<std::string>& unparsed_regions,
                                         const ReferenceGenome& reference)
{
    std::vector<GenomicRegion> result {};
    result.reserve(unparsed_regions.size());
    
    for (const auto& unparsed_region : unparsed_regions) {
        result.push_back(io::parse_region(unparsed_region, reference));
    }
    
    return result;
}

auto transform_to_zero_based(std::vector<GenomicRegion>&& one_based_regions)
{
    std::vector<GenomicRegion> result {};
    result.reserve(one_based_regions.size());
    
    for (auto&& region : one_based_regions) {
        if (region.begin() > 0) {
            result.push_back(shift(std::move(region), -1));
        } else {
            result.push_back(std::move(region));
        }
    }
    
    return result;
}

auto transform_to_zero_based(InputRegionMap::mapped_type&& one_based_regions)
{
    MappableFlatSet<GenomicRegion> result {};
    
    for (auto&& region : one_based_regions) {
        result.insert(shift(std::move(region), -1));
    }
    
    return result;
}

auto transform_to_zero_based(InputRegionMap&& one_based_search_regions)
{
    InputRegionMap result {one_based_search_regions.size()};
    
    for (auto& p : one_based_search_regions) {
        result.emplace(p.first, transform_to_zero_based(std::move(p.second)));
    }
    
    return result;
}

class MissingRegionPathFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_search_regions";
    }
public:
    MissingRegionPathFile(fs::path p) : MissingFileError {std::move(p), "region path"} {};
};

InputRegionMap get_search_regions(const OptionMap& options, const ReferenceGenome& reference)
{
    using namespace utils;
    
    std::vector<GenomicRegion> skip_regions {};
    
    if (options.count("skip-regions") == 1) {
        const auto& region_strings = options.at("skip-regions").as<std::vector<std::string>>();
        append(parse_regions(region_strings, reference), skip_regions);
    }
    
    if (options.count("skip-regions-file") == 1) {
        const auto& input_path = options.at("skip-regions-file").as<std::string>();
        
        auto resolved_path = resolve_path(input_path, options);
        
        if (!fs::exists(resolved_path)) {
            MissingRegionPathFile e {resolved_path};
            e.set_location_specified("the command line option '--skip-regions-file'");
            throw e;
        }
        
        auto regions = io::extract_regions(resolved_path, reference);
        
        if (regions.empty()) {
            logging::WarningLogger log {};
            stream(log) << "The regions path file you specified " << resolved_path
                << " in the command line option '--skip-regions-file' is empty";
        }
        
        append(std::move(regions), skip_regions);
    }
    
    if (options.at("one-based-indexing").as<bool>()) {
        skip_regions = transform_to_zero_based(std::move(skip_regions));
    }
    
    if (options.count("regions") == 0 && options.count("regions-file") == 0) {
        if (options.count("regenotype") == 1) {
            // TODO: only extract regions in the regenotype VCF
        }
        return extract_search_regions(reference, skip_regions);
    }
    
    std::vector<GenomicRegion> input_regions {};
    
    if (options.count("regions") == 1) {
        const auto& region_strings = options.at("regions").as<std::vector<std::string>>();
        append(parse_regions(region_strings, reference), input_regions);
    }
    
    if (options.count("regions-file") == 1) {
        const auto& input_path = options.at("regions-file").as<std::string>();
        
        auto resolved_path = resolve_path(input_path, options);
        
        if (!fs::exists(resolved_path)) {
            MissingRegionPathFile e {resolved_path};
            e.set_location_specified("the command line option '--regions-file'");
            throw e;
        }
        
        auto regions = io::extract_regions(resolved_path, reference);
        
        if (regions.empty()) {
            logging::WarningLogger log {};
            stream(log) << "The regions path file you specified " << resolved_path
                << " in the command line option '--skip-regions-file' is empty";
        }
        
        append(std::move(regions), input_regions);
    }
    
    auto result = extract_search_regions(input_regions, skip_regions);
    
    if (options.at("one-based-indexing").as<bool>()) {
        return transform_to_zero_based(std::move(result));
    }
    
    return result;
}

ContigOutputOrder get_contig_output_order(const OptionMap& options)
{
    return options.at("contig-output-order").as<ContigOutputOrder>();
}

boost::optional<std::vector<SampleName>> get_user_samples(const OptionMap& options)
{
    if (options.count("samples") == 1) {
        return options.at("samples").as<std::vector<SampleName>>();
    }
    return boost::none;
}

class MissingReadPathFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_read_paths";
    }
public:
    MissingReadPathFile(fs::path p) : MissingFileError {std::move(p), "read path"} {};
};

std::vector<fs::path> get_read_paths(const OptionMap& options)
{
    using namespace utils;
    
    std::vector<fs::path> result {};
    
    if (options.count("reads") == 1) {
        auto resolved_paths = resolve_paths(options.at("reads").as<std::vector<std::string>>(), options);
        append(std::move(resolved_paths), result);
    }
    
    if (options.count("reads-file") == 1) {
        const fs::path input_path {options.at("reads-file").as<std::string>()};
        
        auto resolved_path = resolve_path(input_path, options);
        
        if (!fs::exists(resolved_path)) {
            MissingReadPathFile e {resolved_path};
            e.set_location_specified("the command line option '--reads-file'");
            throw e;
        }
        
        auto paths = extract_paths_from_file(resolved_path, options);
        
        auto resolved_paths = resolve_paths(paths, options);
        
        if (resolved_paths.empty()) {
            logging::WarningLogger log {};
            stream(log) << "The read path file you specified " << resolved_path
                        << " in the command line option '--reads-file' is empty";
        }
        
        append(std::move(resolved_paths), result);
    }
    
    std::sort(std::begin(result), std::end(result));
    
    const auto it = std::unique(std::begin(result), std::end(result));
    
    const auto num_duplicates = std::distance(it, std::end(result));
    
    if (num_duplicates > 0) {
        logging::WarningLogger log {};
        stream(log) << "There are " << num_duplicates
                    << " duplicate read paths but only unique paths will be considered";
    }
    
    result.erase(it, std::end(result));
    
    return result;
}

ReadManager make_read_manager(const OptionMap& options)
{
    auto read_paths = get_read_paths(options);
    
    const auto max_open_files = as_unsigned("max-open-read-files", options);
    
    return ReadManager {std::move(read_paths), max_open_files};
}

auto make_read_transformer(const OptionMap& options)
{
    using namespace octopus::readpipe;
    
    ReadTransformer result {};
    
    result.register_transform(CapBaseQualities {125});
    
    if (options.at("disable-read-transforms").as<bool>()) {
        return result;
    }
    
    if (options.count("mask-tails")) {
        const auto tail_mask_size = as_unsigned("mask-tails", options);
        
        if (tail_mask_size > 0) {
            result.register_transform(MaskTail {tail_mask_size});
        }
    }
    
    if (!options.at("disable-soft-clip-masking").as<bool>()) {
        const auto soft_clipped_mask_size = as_unsigned("mask-soft-clipped-boundries", options);
        
        if (soft_clipped_mask_size > 0) {
            result.register_transform(MaskSoftClippedBoundries {soft_clipped_mask_size});
        } else {
            result.register_transform(MaskSoftClipped {});
        }
    }
    
    if (!options.at("disable-adapter-masking").as<bool>()) {
        result.register_transform(MaskAdapters {});
    }
    
    if (!options.at("disable-overlap-masking").as<bool>()) {
        result.register_transform(MaskOverlappedSegment {});
    }
    
    result.shrink_to_fit();
    
    return result;
}

auto make_read_filterer(const OptionMap& options)
{
    using std::make_unique;
    
    using namespace octopus::readpipe;
    
    using ReadFilterer = ReadPipe::ReadFilterer;
    
    ReadFilterer result {};
    
    // these filters are mandatory
    result.add(make_unique<HasValidQualities>());
    result.add(make_unique<HasWellFormedCigar>());
    
    if (options.at("disable-read-filtering").as<bool>()) {
        return result;
    }
    
    if (!options.at("consider-unmapped-reads").as<bool>()) {
        result.add(make_unique<IsMapped>());
    }
    
    const auto min_mapping_quality = as_unsigned("min-mapping-quality", options);
    
    if (min_mapping_quality > 0) {
        result.add(make_unique<IsGoodMappingQuality>(min_mapping_quality));
    }
    
    const auto min_base_quality = as_unsigned("good-base-quality", options);
    const auto min_good_bases   = as_unsigned("min-good-bases", options);
    
    if (min_base_quality > 0 && min_good_bases > 0) {
        result.add(make_unique<HasSufficientGoodQualityBases>(min_base_quality,
                                                                          min_good_bases));
    }
    
    if (min_base_quality > 0 && options.count("min-good-base-fraction") == 1) {
        auto min_good_base_fraction = options.at("min-good-base-fraction").as<double>();
        result.add(make_unique<HasSufficientGoodBaseFraction>(min_base_quality,
                                                                          min_good_base_fraction));
    }
    
    if (options.count("min-read-length") == 1) {
        result.add(make_unique<IsShort>(as_unsigned("min-read-length", options)));
    }
    
    if (options.count("max-read-length") == 1) {
        result.add(make_unique<IsLong>(as_unsigned("max-read-length", options)));
    }
    
    if (!options.at("allow-marked-duplicates").as<bool>()) {
        result.add(make_unique<IsNotMarkedDuplicate>());
    }
    
    if (!options.at("allow-octopus-duplicates").as<bool>()) {
        result.add(make_unique<IsNotDuplicate<ReadFilterer::ReadIterator>>());
    }
    
    if (!options.at("allow-qc-fails").as<bool>()) {
        result.add(make_unique<IsNotMarkedQcFail>());
    }
    
    if (options.at("no-secondary-alignments").as<bool>()) {
        result.add(make_unique<IsNotSecondaryAlignment>());
    }
    
    if (options.at("no-supplementary-alignmenets").as<bool>()) {
        result.add(make_unique<IsNotSupplementaryAlignment>());
    }
    
    if (!options.at("consider-reads-with-unmapped-segments").as<bool>()) {
        result.add(make_unique<IsNextSegmentMapped>());
        result.add(make_unique<IsProperTemplate>());
    }
    
    if (!options.at("consider-reads-with-distant-segments").as<bool>()) {
        result.add(make_unique<IsLocalTemplate>());
    }
    
    if (!options.at("allow-adapter-contaminated-reads").as<bool>()) {
        result.add(make_unique<IsNotContaminated>());
    }
    
    result.shrink_to_fit();
    
    return result;
}

boost::optional<readpipe::Downsampler> make_downsampler(const OptionMap& options)
{
    using namespace octopus::readpipe;
    
    if (options.at("disable-downsampling").as<bool>()) return boost::none;
    
    auto max_coverage    = as_unsigned("downsample-above", options);
    auto target_coverage = as_unsigned("downsample-target", options);
    
    return Downsampler {max_coverage, target_coverage};
}

ReadPipe make_read_pipe(ReadManager& read_manager, std::vector<SampleName> samples,
                        const OptionMap& options)
{
    return ReadPipe {
        read_manager,
        make_read_transformer(options),
        make_read_filterer(options),
        make_downsampler(options),
        std::move(samples)
    };
}

auto make_variant_generator_builder(const OptionMap& options)
{
    logging::WarningLogger warning_log {};
    logging::ErrorLogger log {};
    
    using VGB = coretools::VariantGenerator::Builder;
    
    VGB result {};
    
    if (options.count("generate-candidates-from-source") == 1) {
        result.add_generator(VGB::Generator::External);
        
        const fs::path input_path {options.at("generate-candidates-from-source").as<std::string>()};
        
        auto resolved_path = resolve_path(input_path, options);
        
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--generate-candidates-from-source) does not exist";
        }
        
        result.set_variant_source(std::move(resolved_path));
    }
    
    if (options.count("regenotype") == 1) {
        fs::path regenotype_path {options.at("regenotype").as<std::string>()};
        
        if (options.count("generate-candidates-from-source") == 1) {
            fs::path input_path {options.at("generate-candidates-from-source").as<std::string>()};
            
            if (regenotype_path != input_path) {
                warning_log << "Running in regenotype mode but given a different source variant file";
            } else {
                return result;
            }
        } else {
            result.add_generator(VGB::Generator::External);
        }
        
        auto resolved_path = resolve_path(regenotype_path, options);
        
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << regenotype_path
            << " given in the input option (--generate-candidates-from-source) does not exist";
        }
        
        result.set_variant_source(std::move(resolved_path));
    }
    
    result.set_min_base_quality(as_unsigned("min-base-quality", options));
    result.set_max_variant_size(as_unsigned("max-variant-size", options));
    
    if (options.count("min-supporting-reads")) {
        auto min_supporting_reads = as_unsigned("min-supporting-reads", options);
        
        if (min_supporting_reads == 0) {
            warning_log << "The option --min_supporting_reads was set to 0 - assuming this is a typo and setting to 1";
            ++min_supporting_reads;
        }
        
        result.set_min_supporting_reads(min_supporting_reads);
    } else {
        result.set_min_supporting_reads(2); // TODO: octopus should automatically calculate this
    }
    
    if (!options.at("disable-raw-cigar-candidate-generator").as<bool>()) {
        result.add_generator(VGB::Generator::Alignment);
    }
    
    if (!options.at("disable-assembly-candidate-generator").as<bool>()) {
        result.add_generator(VGB::Generator::Assembler);
        const auto kmer_sizes = options.at("kmer-size").as<std::vector<int>>();
        
        for (const auto k : kmer_sizes) {
            result.add_kmer_size(k);
        }
        
        if (options.count("assembler-mask-base-quality") == 1) {
            result.set_assembler_min_base_quality(as_unsigned("assembler-mask-base-quality", options));
        }
    }
    
    return result;
}

void print_ambiguous_contig_ploidies(const std::vector<ContigPloidy>& contig_ploidies,
                                     const OptionMap& options)
{
    logging::WarningLogger log {};
    
    log << "Ambiguous ploidies found";
    
    for (auto it = std::cbegin(contig_ploidies), end = std::cend(contig_ploidies); it != end;) {
        it = std::adjacent_find(it, std::cend(contig_ploidies),
                                [] (const auto& lhs, const auto& rhs) {
                                    return lhs.contig == rhs.contig;
                                });
        
        if (it != std::cend(contig_ploidies)) {
            const auto it2 = std::find_if(std::next(it), std::cend(contig_ploidies),
                                          [=] (const auto& cp) {
                                              return it->contig != cp.contig;
                                          });
            
            std::ostringstream ss {};
            
            std::copy(it, it2, std::ostream_iterator<ContigPloidy>(ss, " "));
            
            log << ss.str();
            
            it = it2;
        }
    }
}

void remove_duplicate_ploidies(std::vector<ContigPloidy>& contig_ploidies)
{
    std::sort(std::begin(contig_ploidies), std::end(contig_ploidies),
              [] (const auto& lhs, const auto& rhs) {
                  return (lhs.contig == rhs.contig) ? lhs.ploidy < rhs.ploidy : lhs.contig < rhs.contig;
              });
    
    const auto it = std::unique(std::begin(contig_ploidies), std::end(contig_ploidies),
                                [] (const auto& lhs, const auto& rhs) {
                                    return lhs.contig == rhs.contig && lhs.ploidy == rhs.ploidy;
                                });
    
    contig_ploidies.erase(it, std::end(contig_ploidies));
}

bool has_ambiguous_ploidies(const std::vector<ContigPloidy>& contig_ploidies)
{
    const auto it2 = std::adjacent_find(std::cbegin(contig_ploidies), std::cend(contig_ploidies),
                                        [] (const auto& lhs, const auto& rhs) {
                                            return lhs.contig == rhs.contig;
                                        });
    return it2 != std::cend(contig_ploidies);
}

boost::optional<std::vector<ContigPloidy>> extract_contig_ploidies(const OptionMap& options)
{
    std::vector<ContigPloidy> result {};
    
    if (options.count("contig-ploidies-file") == 1) {
        const fs::path input_path {options.at("contig-ploidies-file").as<std::string>()};
        
        const auto resolved_path = resolve_path(input_path, options);
        
        logging::ErrorLogger log {};
        
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--contig-ploidies-file) does not exist";
            return boost::none;
        }
        
        std::ifstream file {resolved_path.string()};
        
        std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                       std::back_inserter(result), [] (const Line& line) {
                           std::istringstream ss {line.line_data};
                           ContigPloidy result {};
                           ss >> result;
                           return result;
                       });
    }
    
    if (options.count("contig-ploidies") == 1) {
        utils::append(options.at("contig-ploidies").as<std::vector<ContigPloidy>>(), result);
    }
    
    remove_duplicate_ploidies(result);
    
    if (has_ambiguous_ploidies(result)) {
        print_ambiguous_contig_ploidies(result, options);
        return boost::none;
    }
    
    return result;
}

bool call_sites_only(const OptionMap& options)
{
    return options.at("sites-only").as<bool>();
}

auto make_haplotype_generator_builder(const OptionMap& options)
{
    using LaggingPolicy = HaplotypeGenerator::Builder::Policies::Lagging;
    
    LaggingPolicy lagging_policy;
    switch (options.at("phasing-level").as<PhasingLevel>()) {
        case PhasingLevel::Minimal:
            lagging_policy = LaggingPolicy::None;
            break;
        case PhasingLevel::Conservative:
            lagging_policy = LaggingPolicy::Conservative;
            break;
        case PhasingLevel::Aggressive:
            lagging_policy = LaggingPolicy::Aggressive;
            break;
    }
    
    const auto max_haplotypes = as_unsigned("max-haplotypes", options);
    
    return HaplotypeGenerator::Builder()
    .set_target_limit(max_haplotypes).set_holdout_limit(2048).set_overflow_limit(16384)
    .set_lagging_policy(lagging_policy).set_max_holdout_depth(3);
}

CallerFactory make_caller_factory(const ReferenceGenome& reference, ReadPipe& read_pipe,
                                  const InputRegionMap& regions, const OptionMap& options)
{
    CallerBuilder vc_builder {
        reference,
        read_pipe,
        make_variant_generator_builder(options),
        make_haplotype_generator_builder(options)
    };
    
    auto caller = options.at("caller").as<std::string>();
    
    if (caller == "population" && read_pipe.num_samples() == 1) {
        caller = "individual";
    }
    
    vc_builder.set_caller(caller);
    
    if (options.count("report-refcalls") == 1) {
        const auto refcall_type = options.at("report-refcalls").as<RefCallType>();
        
        if (refcall_type == RefCallType::Positional) {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::Positional);
        } else {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::Blocked);
        }
    } else {
        vc_builder.set_refcall_type(CallerBuilder::RefCallType::None);
    }
    
    auto min_variant_posterior = options.at("min-variant-posterior").as<Phred<double>>();
    
    if (options.count("regenotype") == 1) {
        if (caller == "cancer") {
            vc_builder.set_min_variant_posterior(min_variant_posterior);
        } else {
            vc_builder.set_min_variant_posterior(Phred<double> {1});
        }
    } else {
        vc_builder.set_min_variant_posterior(min_variant_posterior);
    }
    
    auto min_refcall_posterior = options.at("min-refcall-posterior").as<Phred<double>>();
    
    vc_builder.set_min_refcall_posterior(min_refcall_posterior);
    vc_builder.set_max_haplotypes(as_unsigned("max-haplotypes", options));
    vc_builder.set_min_haplotype_posterior(options.at("min-haplotype-filter-posterior").as<float>());
    
    auto min_phase_score = options.at("min-phase-score").as<Phred<double>>();
    vc_builder.set_min_phase_score(min_phase_score);
    
    vc_builder.set_snp_heterozygosity(options.at("snp-heterozygosity").as<float>());
    vc_builder.set_indel_heterozygosity(options.at("indel-heterozygosity").as<float>());
    
    if (caller == "cancer") {
        if (options.count("normal-sample") == 1) {
            vc_builder.set_normal_sample(options.at("normal-sample").as<std::string>());
        }
        
        vc_builder.set_somatic_mutation_rate(options.at("somatic-mutation-rate").as<float>());
        vc_builder.set_min_somatic_frequency(options.at("min-somatic-frequency").as<float>());
        vc_builder.set_credible_mass(options.at("credible-mass").as<float>());
        
        auto min_somatic_posterior = options.at("min-somatic-posterior").as<Phred<double>>();
        
        vc_builder.set_min_somatic_posterior(min_somatic_posterior);
    } else if (caller == "trio") {
        vc_builder.set_maternal_sample(options.at("maternal-sample").as<std::string>());
        vc_builder.set_paternal_sample(options.at("paternal-sample").as<std::string>());
    }
    
    vc_builder.set_model_filtering(!(options.at("disable-call-filtering").as<bool>()
                                     || options.at("disable-model-filtering").as<bool>()));
    
    const auto contig_ploidies = extract_contig_ploidies(options);
    
    if (!contig_ploidies) {
        // TODO
    }
    
    if (call_sites_only(options)) {
        vc_builder.set_sites_only();
    }
    
    vc_builder.set_flank_scoring(!options.at("disable-inactive-flank-scoring").as<bool>());
    
    CallerFactory result {std::move(vc_builder), as_unsigned("organism-ploidy", options)};
    
    for (const auto& p : regions) {
        const auto it = std::find_if(std::cbegin(*contig_ploidies), std::cend(*contig_ploidies),
                                     [&] (const auto& cp) { return cp.contig == p.first; });
        if (it != std::cend(*contig_ploidies)) {
            result.set_contig_ploidy(p.first, it->ploidy);
        }
    }
    
    return result;
}

fs::path get_final_output_path(const OptionMap& options)
{
    return resolve_path(options.at("output").as<std::string>(), options);
}

VcfWriter make_output_vcf_writer(const OptionMap& options)
{
    return VcfWriter {get_final_output_path(options)};
}

boost::optional<fs::path> create_temp_file_directory(const OptionMap& options)
{
    const auto working_directory = get_working_directory(options);
    
    auto result = working_directory;
    
    const fs::path temp_dir_base_name {"octopus-temp"};
    
    result /= temp_dir_base_name;
    
    constexpr unsigned temp_dir_name_count_limit {10000};
    
    unsigned temp_dir_counter {2};
    
    logging::WarningLogger log {};
    
    while (fs::exists(result) && temp_dir_counter <= temp_dir_name_count_limit) {
        if (fs::is_empty(result)) {
            stream(log) << "Found empty temporary directory " << result
            << ", it may need to be deleted manually";
        }
        result = working_directory;
        result /= temp_dir_base_name.string() + "-" + std::to_string(temp_dir_counter);
        ++temp_dir_counter;
    }
    
    if (temp_dir_counter > temp_dir_name_count_limit) {
        log << "There are many temporary directories in working directory indicating an error"
        " - new directory request blocked";
        return boost::none;
    }
    
    if (!fs::create_directory(result)) {
        stream(log) << "Failed to create temporary directory " << result << " - check permissions";
        return boost::none;
    }
    
    return result;
}

bool legacy_vcf_requested(const OptionMap& options)
{
    return options.at("legacy").as<bool>();
}
} // namespace options
} // namespace octopus
