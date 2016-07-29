//
//  option_collation.cpp
//  Octopus
//
//  Created by Daniel Cooke on 13/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "option_collation.hpp"

#include <string>
#include <iostream>
#include <cctype>
#include <fstream>
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

#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "read_filter.hpp"
#include "read_filterer.hpp"
#include "read_transformer.hpp"
#include "read_transform.hpp"
#include "read_utils.hpp"
#include "candidate_generator_builder.hpp"
#include "haplotype_generator.hpp"
#include "variant_caller_builder.hpp"
#include "variant_caller_factory.hpp"
#include "vcf_reader.hpp"
#include "vcf_writer.hpp"
#include "mappable_algorithms.hpp"
#include "string_utils.hpp"
#include "append.hpp"
#include "phred.hpp"
#include "maths.hpp"
#include "logging.hpp"

namespace octopus { namespace Options
{
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

boost::optional<fs::path> get_home_dir()
{
    static const auto result = fs::path(std::getenv("HOME"));
    
    if (fs::is_directory(result)) {
        return result;
    }
    
    return boost::none;
}

bool is_shorthand_user_path(const fs::path& path)
{
    return !path.empty() && path.string().front() == '~';
}

fs::path expand_user_path(const fs::path& path)
{
    if (is_shorthand_user_path(path)) {
        if (path.string().size() > 1 && path.string()[1] == '/') {
            const auto home_dir = get_home_dir();
            if (home_dir) {
                return fs::path {home_dir->string() + path.string().substr(1)};
            }
            std::ostringstream ss {};
            ss << "Unable to expand user path";
            ss << path;
            ss << " as the user home directory cannot be located";
            throw std::runtime_error {ss.str()};
        }
        return path;
    }
    return path;
}

fs::path get_working_directory(const OptionMap& options)
{
    if (options.count("working-directory") == 1) {
        auto result = expand_user_path(options.at("working-directory").as<std::string>());
        
        if (!fs::exists(result) && !fs::is_directory(result)) {
            std::ostringstream ss {};
            ss << "The working directory ";
            ss << result;
            ss << " given in the option (--working-directory) does not exist";
            throw std::runtime_error {ss.str()};
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

namespace
{
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

std::vector<fs::path>
extract_paths_from_file(const fs::path& file_path, const OptionMap& options)
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
    std::vector<fs::path> good_paths {}, bad_paths {};
    good_paths.reserve(paths.size());
    
    for (const auto& path : paths) {
        try {
            good_paths.push_back(resolve_path(path, options));
        } catch (...) {
            bad_paths.push_back(path);
        }
    }
    
    return std::make_pair(std::move(good_paths), std::move(bad_paths));
}

auto resolve_paths(const std::vector<std::string>& path_strings,
                   const OptionMap& options)
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
        num_threads = options.at("threads").as<unsigned>();
    }
    
    return num_threads != 1;
}

boost::optional<unsigned> get_num_threads(const OptionMap& options)
{
    unsigned num_threads {1};
    
    if (options.count("threads") == 1) {
        num_threads = options.at("threads").as<unsigned>();
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
    Logging::ErrorLogger log {};
    
    const fs::path input_path {options.at("reference").as<std::string>()};
    
    auto resolved_path = resolve_path(input_path, options);
    
    if (!fs::exists(resolved_path)) {
        stream(log) << "The path " << input_path
        << " given in the input option (--reference) does not exist";
    }
    
    if (!is_file_readable(resolved_path)) {
        stream(log) << "The path " << input_path
        << " given in the input option (--reference) is not readable";
    }
    
    const auto ref_cache_size = options.at("max-reference-cache-footprint").as<float>();
    
    static constexpr unsigned Scale {1'000'000};
    
    return ::make_reference(std::move(resolved_path),
                            static_cast<std::size_t>(Scale * ref_cache_size),
                            is_threading_allowed(options));
}

bool is_bed_file(const fs::path& path)
{
    return path.extension().string() == ".bed";
}

void seek_past_bed_header(std::ifstream& bed_file)
{
    // TODO
}

std::string convert_bed_line_to_region_str(const std::string& bed_line)
{
    constexpr static char bed_delim {'\t'};
    
    const auto tokens = split(bed_line, bed_delim);
    
    switch (tokens.size()) {
        case 0:
            throw std::runtime_error {"BadBED: found empty BED record"};
        case 1:
            return std::string {tokens[0]};
        case 2:
            // Assume this represents a half range rather than a position
            return std::string {tokens[0] + ':' + tokens[1] + '-'};
        default:
            return std::string {tokens[0] + ':' + tokens[1] + '-' + tokens[2]};
    }
}

std::function<GenomicRegion(const std::string&)>
make_region_line_parser(const fs::path& region_path, const ReferenceGenome& reference)
{
    if (is_bed_file(region_path)) {
        return [&] (const std::string& line) -> GenomicRegion
        {
            return parse_region(convert_bed_line_to_region_str(line), reference);
        };
    } else {
        return [&] (const std::string& line) { return parse_region(line, reference); };
    }
}

auto extract_regions_from_file(const fs::path& file_path, const ReferenceGenome& reference)
{
    std::ifstream file {file_path.string()};
    
    if (is_bed_file(file_path)) {
        seek_past_bed_header(file);
    }
    
    std::deque<GenomicRegion> result {};
    
    std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                   std::back_inserter(result), make_region_line_parser(file_path, reference));
    
    result.shrink_to_fit();
    
    return result;
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

MappableFlatSet<GenomicRegion>
cut(const MappableFlatSet<GenomicRegion>& mappables, const MappableFlatSet<GenomicRegion>& regions)
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
    
    bool all_region_parsed {true};
    
    for (const auto& unparsed_region : unparsed_regions) {
        Logging::WarningLogger log {};
        try {
            result.push_back(parse_region(unparsed_region, reference));
        } catch (std::exception& e) {
            all_region_parsed = false;
            stream(log) << "Could not parse input region \"" << unparsed_region
            << "\". Check the format is correct, the contig is in the reference genome \""
            << reference.name() << "\", and the coordinate range is in bounds.";
        }
    }
    
    if (!all_region_parsed) {
        result.clear();
        result.shrink_to_fit();
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

InputRegionMap get_search_regions(const OptionMap& options, const ReferenceGenome& reference)
{
    Logging::ErrorLogger log {};
    
    std::vector<GenomicRegion> skip_regions {};
    
    bool all_parsed {true};
    
    if (options.count("skip-regions") == 1) {
        const auto& region_strings = options.at("skip-regions").as<std::vector<std::string>>();
        auto parsed_regions = parse_regions(region_strings, reference);
        if (region_strings.size() == parsed_regions.size()) {
            append(std::move(parsed_regions), skip_regions);
        } else {
            all_parsed = false;
        }
    }
    
    if (options.count("skip-regions-file") == 1) {
        const auto& input_path = options.at("skip-regions-file").as<std::string>();
        
        auto resolved_path = resolve_path(input_path, options);
        
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--skip-regions-file) does not exist";
        } else if (!is_file_readable(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--skip-regions-file) is not readable";
        } else {
            append(extract_regions_from_file(resolved_path, reference), skip_regions);
        }
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
        auto parsed_regions = parse_regions(region_strings, reference);
        if (region_strings.size() == parsed_regions.size()) {
            append(std::move(parsed_regions), input_regions);
        } else {
            all_parsed = false;
        }
    }
    
    if (options.count("regions-file") == 1) {
        const auto& input_path = options.at("regions-file").as<std::string>();
        
        auto resolved_path = resolve_path(input_path, options);
        
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--skip-regions-file) does not exist";
        } else if (!is_file_readable(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--skip-regions-file) is not readable";
        } else {
            append(extract_regions_from_file(resolved_path, reference), input_regions);
        }
    }
    
    if (!all_parsed) {
        Logging::WarningLogger log {};
        if (!input_regions.empty()) {
            stream(log) << "Detected unparsed input regions so dumping "
            << input_regions.size() << " parsed regions";
            input_regions.clear();
        }
        skip_regions.clear();
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

namespace
{
    void log_unresolved_read_paths(const std::vector<fs::path>& paths,
                                   const std::string& option)
    {
        Logging::WarningLogger log {};
        for (const auto& path : paths) {
            stream(log) << "Could not resolve the path " << path
            << " given in the input option (--" + option +")";
        }
    }
    
    auto parition_existent_paths(std::vector<fs::path>& paths)
    {
        return std::partition(std::begin(paths), std::end(paths),
                              [] (const auto& path) { return fs::exists(path); });
    }
    
    template <typename InputIt>
    void log_nonexistent_read_paths(InputIt first, InputIt last, const std::string& option)
    {
        Logging::WarningLogger log {};
        std::for_each(first, last, [&option, &log] (const auto& path) {
            stream(log) << "The path " << path
            << " given in the input option (--" + option + ") does not exist";
        });
    }
    
    auto parition_readable_paths(std::vector<fs::path>& paths)
    {
        return std::partition(std::begin(paths), std::end(paths),
                              [] (const auto& path) { return is_file_readable(path); });
    }
    
    template <typename InputIt>
    void log_unreadable_read_paths(InputIt first, InputIt last, const std::string& option)
    {
        Logging::WarningLogger log {};
        std::for_each(first, last, [&option, &log] (const auto& path) {
            stream(log) << "The path " << path
            << " given in the input option (--" + option + ") is not readable";
        });
    }
} // namespace

boost::optional<std::vector<fs::path>> get_read_paths(const OptionMap& options)
{
    Logging::ErrorLogger log {};
    
    std::vector<fs::path> result {};
    
    bool all_paths_good {true};
    
    std::vector<fs::path> resolved_paths {}, unresolved_paths {};
    
    if (options.count("reads") == 1) {
        const auto& read_paths = options.at("reads").as<std::vector<std::string>>();
        
        std::tie(resolved_paths, unresolved_paths) = resolve_paths(read_paths, options);
        
        if (!unresolved_paths.empty()) {
            log_unresolved_read_paths(unresolved_paths, "reads");
            all_paths_good = false;
        }
        
        auto it = parition_existent_paths(resolved_paths);
        
        if (it != std::end(resolved_paths)) {
            log_nonexistent_read_paths(it, std::end(resolved_paths), "reads");
            all_paths_good = false;
            resolved_paths.erase(it, std::end(resolved_paths));
        }
        
        it = parition_readable_paths(resolved_paths);
        
        if (it != std::end(resolved_paths)) {
            log_unreadable_read_paths(it, std::end(resolved_paths), "reads");
            all_paths_good = false;
        }
        
        append(std::move(resolved_paths), result);
    }
    
    if (options.count("reads-file") == 1) {
        // first we need to make sure the path to the paths is okay
        const fs::path input_path {options.at("reads-file").as<std::string>()};
        
        auto resolved_path = resolve_path(input_path, options);
        
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--reads-file) does not exist";
            all_paths_good = false;
        } else if (!is_file_readable(resolved_path)) {
            stream(log) << "The path " << input_path
            << " given in the input option (--reads-file) is not readable";
            all_paths_good = false;
        } else {
            auto paths = extract_paths_from_file(resolved_path, options);
            
            std::tie(resolved_paths, unresolved_paths) = resolve_paths(paths, options);
            
            if (!unresolved_paths.empty()) {
                log_unresolved_read_paths(unresolved_paths, "reads-file");
                all_paths_good = false;
            }
            
            auto it = parition_existent_paths(resolved_paths);
            
            if (it != std::end(resolved_paths)) {
                log_nonexistent_read_paths(it, std::end(resolved_paths), "reads-file");
                all_paths_good = false;
                resolved_paths.erase(it, std::end(resolved_paths));
            }
            
            it = parition_readable_paths(resolved_paths);
            
            if (it != std::end(resolved_paths)) {
                log_unreadable_read_paths(it, std::end(resolved_paths), "reads-file");
                all_paths_good = false;
            }
            
            append(std::move(resolved_paths), result);
        }
    }
    
    std::sort(std::begin(result), std::end(result));
    
    const auto it = std::unique(std::begin(result), std::end(result));
    
    const auto num_duplicates = std::distance(it, std::end(result));
    
    if (num_duplicates > 0) {
        Logging::WarningLogger log {};
        stream(log) << "There are " << num_duplicates
        << " duplicate read paths but only unique paths will be considered";
    }
    
    result.erase(it, std::end(result));
    
    if (!all_paths_good && result.size() > 0) {
        Logging::WarningLogger log {};
        auto slog = stream(log);
        slog << "There are bad read paths so dumping " << result.size() << " good path";
        if (result.size() > 1) slog << "s";
        result.clear();
    }
    
    return result;
}

ReadManager make_read_manager(const OptionMap& options)
{
    auto read_paths = get_read_paths(options);
    
    if (read_paths) {
        const auto max_open_files = options.at("max-open-read-files").as<unsigned>();
        
        return ReadManager {*std::move(read_paths), max_open_files};
    }
    
    throw std::runtime_error {"Unable to load read paths"};
}

ReadTransformer make_read_transformer(const OptionMap& options)
{
    using namespace read_transform;
    
    ReadTransformer result {};
    
    result.register_transform(CapBaseQualities {125});
    
    if (options.at("disable-read-transforms").as<bool>()) {
        return result;
    }
    
    if (options.count("mask-tails")) {
        const auto tail_mask_size = options.at("mask-tails").as<unsigned>();
        
        if (tail_mask_size > 0) {
            result.register_transform(MaskTail {tail_mask_size});
        }
    }
    
    if (!options.at("disable-soft-clip-masking").as<bool>()) {
        const auto soft_clipped_mask_size = options.at("mask-soft-clipped-boundries").as<unsigned>();
        
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

ReadPipe::ReadFilterer make_read_filterer(const OptionMap& options)
{
    using std::make_unique;
    
    using namespace read_filter;
    
    using ReadFilterer = ReadPipe::ReadFilterer;
    
    ReadFilterer result {};
    
    // these filters are mandatory
    result.register_filter(make_unique<HasValidQualities>());
    result.register_filter(make_unique<HasWellFormedCigar>());
    
    if (options.at("disable-read-filtering").as<bool>()) {
        return result;
    }
    
    if (!options.at("consider-unmapped-reads").as<bool>()) {
        result.register_filter(make_unique<IsMapped>());
    }
    
    const auto min_mapping_quality = options.at("min-mapping-quality").as<unsigned>();
    
    if (min_mapping_quality > 0) {
        result.register_filter(make_unique<IsGoodMappingQuality>(min_mapping_quality));
    }
    
    const auto min_base_quality = options.at("good-base-quality").as<unsigned>();
    const auto min_good_bases   = options.at("min-good-bases").as<unsigned>();
    
    if (min_base_quality > 0 && min_good_bases > 0) {
        result.register_filter(make_unique<HasSufficientGoodQualityBases>(min_base_quality,
                                                                          min_good_bases));
    }
    
    if (min_base_quality > 0 && options.count("min-good-base-fraction") == 1) {
        auto min_good_base_fraction = options.at("min-good-base-fraction").as<double>();
        result.register_filter(make_unique<HasSufficientGoodBaseFraction>(min_base_quality,
                                                                          min_good_base_fraction));
    }
    
    if (options.count("min-read-length") == 1) {
        result.register_filter(make_unique<IsShort>(options.at("min-read-length").as<unsigned>()));
    }
    
    if (options.count("max-read-length") == 1) {
        result.register_filter(make_unique<IsLong>(options.at("max-read-length").as<unsigned>()));
    }
    
    if (!options.at("allow-marked-duplicates").as<bool>()) {
        result.register_filter(make_unique<IsNotMarkedDuplicate>());
    }
    
    if (!options.at("allow-octopus-duplicates").as<bool>()) {
        result.register_filter(make_unique<IsNotDuplicate<ReadFilterer::ReadIterator>>());
    }
    
    if (!options.at("allow-qc-fails").as<bool>()) {
        result.register_filter(make_unique<IsNotMarkedQcFail>());
    }
    
    if (options.at("no-secondary-alignments").as<bool>()) {
        result.register_filter(make_unique<IsNotSecondaryAlignment>());
    }
    
    if (options.at("no-supplementary-alignmenets").as<bool>()) {
        result.register_filter(make_unique<IsNotSupplementaryAlignment>());
    }
    
    if (!options.at("consider-reads-with-unmapped-segments").as<bool>()) {
        result.register_filter(make_unique<IsNextSegmentMapped>());
        result.register_filter(make_unique<IsProperTemplate>());
    }
    
    if (!options.at("consider-reads-with-distant-segments").as<bool>()) {
        result.register_filter(make_unique<IsLocalTemplate>());
    }
    
    if (!options.at("allow-adapter-contaminated-reads").as<bool>()) {
        result.register_filter(make_unique<IsNotContaminated>());
    }
    
    result.shrink_to_fit();
    
    return result;
}

boost::optional<Downsampler> make_downsampler(const OptionMap& options)
{
    if (options.at("disable-downsampling").as<bool>()) {
        return boost::none;
    }
    
    auto max_coverage    = options.at("downsample-above").as<unsigned>();
    auto target_coverage = options.at("downsample-target").as<unsigned>();
    
    return Downsampler(max_coverage, target_coverage);
}

CandidateGeneratorBuilder make_candidate_generator_builder(const OptionMap& options,
                                                           const ReferenceGenome& reference)
{
    Logging::WarningLogger warning_log {};
    Logging::ErrorLogger log {};
    
    CandidateGeneratorBuilder result {};
    
    result.set_reference(reference);
    
    if (options.count("generate-candidates-from-source") == 1) {
        result.add_generator(CandidateGeneratorBuilder::Generator::External);
        
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
            result.add_generator(CandidateGeneratorBuilder::Generator::External);
        }
        
        auto resolved_path = resolve_path(regenotype_path, options);
        
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << regenotype_path
            << " given in the input option (--generate-candidates-from-source) does not exist";
        }
        
        result.set_variant_source(std::move(resolved_path));
    }
    
    result.set_min_base_quality(options.at("min-base-quality").as<unsigned>());
    result.set_max_variant_size(options.at("max-variant-size").as<unsigned>());
    
    if (options.count("min-supporting-reads")) {
        auto min_supporting_reads = options.at("min-supporting-reads").as<unsigned>();
        
        if (min_supporting_reads == 0) {
            warning_log << "The option --min_supporting_reads was set to 0 - assuming this is a typo and setting to 1";
            ++min_supporting_reads;
        }
        
        result.set_min_supporting_reads(min_supporting_reads);
    } else {
        result.set_min_supporting_reads(2); // TODO: Octopus should automatically calculate this
    }
    
    if (!options.at("disable-raw-cigar-candidate-generator").as<bool>()) {
        result.add_generator(CandidateGeneratorBuilder::Generator::Alignment);
    }
    
    if (!options.at("disable-assembly-candidate-generator").as<bool>()) {
        result.add_generator(CandidateGeneratorBuilder::Generator::Assembler);
        const auto kmer_sizes = options.at("kmer-size").as<std::vector<unsigned>>();
        
        for (const auto k : kmer_sizes) {
            result.add_kmer_size(k);
        }
        
        if (options.count("assembler-mask-base-quality") == 1) {
            result.set_assembler_min_base_quality(options.at("assembler-mask-base-quality").as<unsigned>());
        }
    }
    
    return result;
}

void print_ambiguous_contig_ploidies(const std::vector<ContigPloidy>& contig_ploidies,
                                     const OptionMap& options)
{
    Logging::WarningLogger log {};
    
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
        
        Logging::ErrorLogger log {};
        
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
        append(options.at("contig-ploidies").as<std::vector<ContigPloidy>>(), result);
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
    
    const auto max_haplotypes = options.at("max-haplotypes").as<unsigned>();
    
    return HaplotypeGenerator::Builder()
    .set_target_limit(max_haplotypes).set_holdout_limit(2048).set_overflow_limit(16384)
    .set_lagging_policy(lagging_policy).set_max_holdout_depth(3);
}

VariantCallerFactory
make_variant_caller_factory(const ReferenceGenome& reference,
                            ReadPipe& read_pipe,
                            const CandidateGeneratorBuilder& candidate_generator_builder,
                            const InputRegionMap& regions,
                            const OptionMap& options)
{
    VariantCallerBuilder vc_builder {
        reference, read_pipe, candidate_generator_builder,
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
            vc_builder.set_refcall_type(VariantCallerBuilder::RefCallType::Positional);
        } else {
            vc_builder.set_refcall_type(VariantCallerBuilder::RefCallType::Blocked);
        }
    } else {
        vc_builder.set_refcall_type(VariantCallerBuilder::RefCallType::None);
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
    vc_builder.set_max_haplotypes(options.at("max-haplotypes").as<unsigned>());
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
    
    VariantCallerFactory result {std::move(vc_builder), options.at("organism-ploidy").as<unsigned>()};
    
    for (const auto& p : regions) {
        const auto it = std::find_if(std::cbegin(*contig_ploidies), std::cend(*contig_ploidies),
                                     [&] (const auto& cp) { return cp.contig == p.first; });
        if (it != std::cend(*contig_ploidies)) {
            result.set_contig_ploidy(p.first, it->ploidy);
        }
    }
    
    return result;
}

boost::optional<fs::path> get_final_output_path(const OptionMap& options)
{
    Logging::ErrorLogger log {};
    
    const auto input_path = options.at("output").as<std::string>();
    
    if (input_path == "-") {
        return fs::path {input_path}; // Output goes to stdout
    }
    
    const auto resolved_path = resolve_path(input_path, options);
    
    if (!is_file_writable(resolved_path)) {
        stream(log) << "The path " << input_path << " given in the input option output is not writable";
        return boost::none;
    }
    
    return resolved_path;
}

VcfWriter make_output_vcf_writer(const OptionMap& options)
{
    auto out_path = get_final_output_path(options);
    
    if (out_path) {
        return VcfWriter {*std::move(out_path)};
    }
    
    return VcfWriter {};
}

boost::optional<fs::path> create_temp_file_directory(const OptionMap& options)
{
    const auto working_directory = get_working_directory(options);
    
    auto result = working_directory;
    
    const fs::path temp_dir_base_name {"octopus-temp"};
    
    result /= temp_dir_base_name;
    
    constexpr unsigned temp_dir_name_count_limit {10000};
    
    unsigned temp_dir_counter {2};
    
    Logging::WarningLogger log {};
    
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
} // namespace Options
} // namespace octopus
