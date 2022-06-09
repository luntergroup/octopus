// Copyright (c) 2015-2021 Daniel Cooke
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

#include "utils/path_utils.hpp"
#include "utils/read_stats.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/string_utils.hpp"
#include "utils/repeat_finder.hpp"
#include "utils/append.hpp"
#include "utils/maths.hpp"
#include "basics/phred.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "basics/trio.hpp"
#include "basics/pedigree.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "core/tools/coretools.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "core/models/error/error_model_factory.hpp"
#include "core/callers/caller_builder.hpp"
#include "logging/logging.hpp"
#include "io/region/region_parser.hpp"
#include "io/pedigree/pedigree_reader.hpp"
#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"
#include "exceptions/user_error.hpp"
#include "exceptions/program_error.hpp"
#include "exceptions/system_error.hpp"
#include "exceptions/missing_file_error.hpp"
#include "exceptions/malformed_file_error.hpp"
#include "core/csr/filters/threshold_filter_factory.hpp"
#include "core/csr/filters/training_filter_factory.hpp"
#include "core/csr/filters/random_forest_filter_factory.hpp"

namespace octopus { namespace options {

bool is_set(const std::string& option, const OptionMap& options) noexcept
{
    return options.count(option) == 1;
}

boost::optional<fs::path> get_output_path(const OptionMap& options);

// unsigned are banned from the option map to prevent user input errors, but once the option
// map is passed they are all safe
unsigned as_unsigned(const std::string& option, const OptionMap& options)
{
    return static_cast<unsigned>(options.at(option).as<int>());
}

bool is_run_command(const OptionMap& options)
{
    return !is_set("help", options) && !is_set("version", options);
}

bool is_debug_mode(const OptionMap& options)
{
    return is_set("debug", options);
}

bool is_trace_mode(const OptionMap& options)
{
    return is_set("trace", options);
}

void emit_in_development_warning(const std::string& option)
{
    logging::WarningLogger log {};
    stream(log) << "The requested option '--" << option 
                << "' invokes a feature that is currently under development"
                    " and may not function correctly or as expected";
}

namespace {

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
    if (is_set("working-directory", options)) {
        auto result = expand_user_path(options.at("working-directory").as<fs::path>());
        if (!fs::exists(result) && !fs::is_directory(result)) {
            throw InvalidWorkingDirectory {result};
        }
        return result;
    }
    return fs::current_path();
}

SymblinkResolvePolicy get_symlink_resolve_policy(const OptionMap& options)
{
    return options.at("resolve-symlinks").as<bool>() ? SymblinkResolvePolicy::resolve : SymblinkResolvePolicy::dont_resolve;
}

fs::path resolve_path(const fs::path& path, const OptionMap& options)
{
    const auto wd_policy = WorkingDirectoryResolvePolicy::prefer_working_directory;
    const auto symlink_policy = get_symlink_resolve_policy(options);
    return ::octopus::resolve_path(path, get_working_directory(options), wd_policy, symlink_policy);
}

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

std::vector<fs::path> extract_paths_from_file(const fs::path& file_path)
{
    std::ifstream file {file_path.string()};
    assert(file.good());
    std::vector<fs::path> result {};
    std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                   std::back_inserter(result), [] (const auto& line) { return line.line_data; });
    result.erase(std::remove_if(std::begin(result), std::end(result), [] (const auto& path) { return path.empty(); }),
                 std::end(result));
    return result;
}

auto resolve_file_paths(const fs::path& file_path, std::vector<fs::path> paths_in_file, const OptionMap& options)
{
    for (auto& path : paths_in_file) {
        if (!fs::exists(path)) {
            auto full_path = file_path.parent_path();
            full_path /= path;
            if (fs::exists(full_path)) {
                path = std::move(full_path);
            } else {
                path = resolve_path(path, options);
            }
        }
    }
    return paths_in_file;
}

auto get_resolved_paths_from_file(const fs::path& file, const OptionMap& options)
{
    return resolve_file_paths(file, extract_paths_from_file(file), options);
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
    const bool result {test.is_open()};
    fs::remove(path);
    return result;
}

} // namespace

bool is_threading_allowed(const OptionMap& options)
{
    unsigned num_threads {1};
    if (is_set("threads", options)) {
        num_threads = as_unsigned("threads", options);
    }
    return num_threads != 1;
}

boost::optional<unsigned> get_num_threads(const OptionMap& options)
{
    unsigned num_threads {1};
    if (is_set("threads", options)) {
        num_threads = as_unsigned("threads", options);
    }
    if (num_threads > 0) return num_threads;
    return boost::none;
}

ExecutionPolicy get_thread_execution_policy(const OptionMap& options)
{
    if (is_set("threads", options)) {
        if (options.at("threads").as<int>() == 0) {
            return ExecutionPolicy::par;
        } else {
            return ExecutionPolicy::seq;
        }
    } else {
        return ExecutionPolicy::seq;
    }
}

MemoryFootprint get_target_read_buffer_size(const OptionMap& options)
{
    return options.at("target-read-buffer-memory").as<MemoryFootprint>();
}

boost::optional<fs::path> get_debug_log_file_name(const OptionMap& options)
{
    if (is_debug_mode(options)) {
        return resolve_path(options.at("debug").as<fs::path>(), options);
    } else {
        return boost::none;
    }
}

boost::optional<fs::path> get_trace_log_file_name(const OptionMap& options)
{
    if (is_trace_mode(options)) {
        return resolve_path(options.at("trace").as<fs::path>(), options);
    } else {
        return boost::none;
    }
}

bool is_fast_mode(const OptionMap& options)
{
    return options.at("fast").as<bool>() || options.at("very-fast").as<bool>();
}

bool is_very_fast_mode(const OptionMap& options)
{
    return options.at("very-fast").as<bool>();
}

ReferenceGenome make_reference(const OptionMap& options)
{
    const fs::path input_path {options.at("reference").as<fs::path>()};
    auto resolved_path = resolve_path(input_path, options);
    auto ref_cache_size = options.at("max-reference-cache-memory").as<MemoryFootprint>();
    static constexpr MemoryFootprint min_non_zero_reference_cache_size {1'000}; // 1Kb
    if (ref_cache_size.bytes() > 0 && ref_cache_size < min_non_zero_reference_cache_size) {
        static bool warned {false};
        if (!warned) {
            logging::WarningLogger warn_log {};
            stream(warn_log) << "Ignoring given reference cache size of " << ref_cache_size
                             << " as this size is too small. The maximum cache size will be set to zero";
            warned = true;
        }
        ref_cache_size = 0;
    }
    static constexpr MemoryFootprint min_warn_non_zero_reference_cache_size {1'000'000}; // 1Mb
    if (ref_cache_size.bytes() > 0 && ref_cache_size < min_warn_non_zero_reference_cache_size) {
        static bool warned {false};
        if (!warned) {
            logging::WarningLogger warn_log {};
            stream(warn_log) << "The given reference cache size " << ref_cache_size
                             << " is very small and may not result in good performance.";
            warned = true;
        }
    }
    try {
        return octopus::make_reference(std::move(resolved_path), ref_cache_size, is_threading_allowed(options));
    } catch (MissingFileError& e) {
        e.set_location_specified("the command line option --reference");
        throw;
    } catch (...) {
        throw;
    }
}

InputRegionMap make_search_regions(const std::vector<GenomicRegion>& regions)
{
    std::map<ContigName, std::deque<GenomicRegion>> contig_mapped_regions {};
    for (const auto& region : regions) {
        contig_mapped_regions[region.contig_name()].push_back(region);
    }
    InputRegionMap result {};
    result.reserve(contig_mapped_regions.size());
    for (auto& p : contig_mapped_regions) {
        std::sort(std::begin(p.second), std::end(p.second));
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

auto get_unskipped(const MappableFlatSet<GenomicRegion>& regions, const MappableFlatSet<GenomicRegion>& skips)
{
    if (skips.empty()) return regions;
    MappableFlatSet<GenomicRegion> result {};
    for (const auto& region : regions) {
        const auto overlapped = skips.overlap_range(region);
        if (empty(overlapped)) {
            result.emplace(region);
        } else if (!contains(overlapped.front(), region)) {
            if (begins_before(region, overlapped.front())) {
                result.emplace(left_overhang_region(region, overlapped.front()));
            }
            auto intervening_chunks = extract_intervening_regions(overlapped);
            using std::make_move_iterator;
            result.insert(make_move_iterator(std::begin(intervening_chunks)),
                          make_move_iterator(std::end(intervening_chunks)));
            if (ends_before(overlapped.back(), region)) {
                result.emplace(right_overhang_region(region, overlapped.back()));
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
            result.emplace(p.first, get_unskipped(std::move(p.second), skipped.at(p.first)));
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

bool is_region_range(const std::vector<std::string>& unparsed_regions, const ReferenceGenome& reference)
{
    return unparsed_regions.size() == 3 && unparsed_regions[1] == "to" && !reference.has_contig("to");
}

class BadRegionRange : public UserError
{
    std::string do_where() const override
    {
        return "make_region_range";
    }
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The region " << lhs_ << " is after " << rhs_;
        if (!is_same_contig(lhs_, rhs_)) {
            ss << " in reference index";
        }
        return ss.str();
    }
    std::string do_help() const override
    {
        return "Ensure the region range format is <lhs> to <rhs> where lhs occurs before rhs in the reference index";
    }
    
    GenomicRegion lhs_, rhs_;
public:
    BadRegionRange(GenomicRegion lhs, GenomicRegion rhs) : lhs_ {lhs}, rhs_ {rhs} {}
};

std::vector<GenomicRegion> make_region_range(GenomicRegion lhs, GenomicRegion rhs, const ReferenceGenome& reference)
{
    assert(reference.has_contig(lhs.contig_name()) && reference.has_contig(rhs.contig_name()));
    std::vector<GenomicRegion> result {};
    if (is_same_contig(lhs, rhs)) {
        if (lhs == rhs || (begins_before(lhs, rhs) && ends_before(lhs, rhs))) {
            result.push_back(closed_region(lhs, rhs));
        } else {
            throw BadRegionRange {lhs, rhs};
        }
    } else {
        auto reference_contigs = reference.contig_names();
        const auto lhs_itr = std::find(std::cbegin(reference_contigs), std::cend(reference_contigs), lhs.contig_name());
        assert(lhs_itr != std::cend(reference_contigs));
        const auto rhs_itr = std::find(std::next(lhs_itr), std::cend(reference_contigs), rhs.contig_name());
        if (rhs_itr != std::cend(reference_contigs)) {
            result.reserve(std::distance(lhs_itr, rhs_itr) + 1);
            result.emplace_back(lhs.contig_name(), lhs.begin(), reference.contig_region(lhs.contig_name()).end());
            std::transform(std::next(lhs_itr), rhs_itr, std::back_inserter(result),
                           [&] (const auto& contig) { return reference.contig_region(contig); });
            result.emplace_back(rhs.contig_name(), 0, rhs.end());
        } else {
            throw BadRegionRange {lhs, rhs};
        }
    }
    return result;
}

std::vector<GenomicRegion> parse_region_range(const std::string& lhs, const std::string& rhs, const ReferenceGenome& reference)
{
    return make_region_range(io::parse_region(lhs, reference), io::parse_region(rhs, reference), reference);
}

std::vector<GenomicRegion> parse_regions(const std::vector<std::string>& unparsed_regions, const ReferenceGenome& reference)
{
    std::vector<GenomicRegion> result {};
    if (is_region_range(unparsed_regions, reference)) {
        result = parse_region_range(unparsed_regions.front(), unparsed_regions.back(), reference);
    } else {
        result.reserve(unparsed_regions.size());
        for (const auto& unparsed_region : unparsed_regions) {
            result.push_back(io::parse_region(unparsed_region, reference));
        }
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
    if (is_set("skip-regions", options)) {
        const auto& region_strings = options.at("skip-regions").as<std::vector<std::string>>();
        append(parse_regions(region_strings, reference), skip_regions);
    }
    if (is_set("skip-regions-file", options)) {
        const auto& input_path = options.at("skip-regions-file").as<fs::path>();
        auto resolved_path = resolve_path(input_path, options);
        if (!fs::exists(resolved_path)) {
            MissingRegionPathFile e {resolved_path};
            e.set_location_specified("the command line option '--skip-regions-file'");
            throw e;
        }
        auto regions = io::extract_regions(resolved_path, reference, io::NonreferenceContigPolicy::ignore);
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
    if (!is_set("regions", options) && !is_set("regions-file", options)) {
        if (is_set("regenotype", options)) {
            // TODO: only extract regions in the regenotype VCF
        }
        return extract_search_regions(reference, skip_regions);
    }
    std::vector<GenomicRegion> input_regions {};
    if (is_set("regions", options)) {
        const auto& region_strings = options.at("regions").as<std::vector<std::string>>();
        append(parse_regions(region_strings, reference), input_regions);
    }
    if (is_set("regions-file", options)) {
        const auto& input_path = options.at("regions-file").as<fs::path>();
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

bool ignore_unmapped_contigs(const OptionMap& options)
{
    return options.at("ignore-unmapped-contigs").as<bool>();
}

class MissingSamplesFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_user_samples";
    }
public:
    MissingSamplesFile(fs::path p) : MissingFileError {std::move(p), "samples file"} {};
};

std::vector<SampleName> read_user_samples(const fs::path& samples_filename)
{
    if (!fs::exists(samples_filename)) {
        MissingSamplesFile e {samples_filename};
        e.set_location_specified("the command line option '--samples-file'");
        throw e;
    }
    std::ifstream samples_file {samples_filename.string()};
    std::vector<SampleName> result {};
    std::string sample {};
    while (std::getline(samples_file, sample)) {
        if (!sample.empty()) result.push_back(std::move(sample));
    }
    return result;
}

boost::optional<std::vector<SampleName>> get_user_samples(const OptionMap& options)
{
    boost::optional<std::vector<SampleName>> result {};
    if (is_set("samples", options)) {
        result = options.at("samples").as<std::vector<SampleName>>();
    }
    if (is_set("samples-file", options)) {
        if (!result) result = std::vector<SampleName> {};
        const auto samples_filename = resolve_path(options.at("samples-file").as<fs::path>(), options);
        utils::append(read_user_samples(samples_filename), *result);
    }
    if (result) {
        std::sort(std::begin(*result), std::end(*result));
        result->erase(std::unique(std::begin(*result), std::end(*result)), std::end(*result));
    }
    return result;
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

class MalfordReadPathFile : public MalformedFileError
{
    std::string do_where() const override
    {
        return "get_read_paths";
    }
public:
    MalfordReadPathFile(fs::path p) : MalformedFileError {std::move(p), "text"} {};
};

void remove_duplicates(std::vector<fs::path>& paths, const std::string& type, const bool log = true)
{
    std::sort(std::begin(paths), std::end(paths));
    const auto first_duplicate = std::adjacent_find(std::begin(paths), std::end(paths));
    if (first_duplicate != std::end(paths)) {
        std::deque<fs::path> duplicates {};
        for (auto duplicate_itr = first_duplicate; duplicate_itr != std::end(paths); ) {
            duplicates.push_back(*duplicate_itr);
            duplicate_itr = std::adjacent_find(std::find_if(std::next(duplicate_itr, 2), std::end(paths),
                                                            [=] (const auto& path) { return path != *duplicate_itr; }),
                                               std::end(paths));
        }
        const auto num_paths = paths.size();
        paths.erase(std::unique(first_duplicate, std::end(paths)), std::end(paths));
        const auto num_unique_paths = paths.size();
        const auto num_duplicate_paths = num_paths - num_unique_paths;
        if (log) {
            logging::WarningLogger warn_log {};
            auto warn_log_stream = stream(warn_log);
            warn_log_stream << "Ignoring " << num_duplicate_paths << " duplicate " << type << " path";
            if (num_duplicate_paths > 1) {
                warn_log_stream << 's';
            }
            warn_log_stream << ": ";
            std::for_each(std::cbegin(duplicates), std::prev(std::cend(duplicates)), [&] (const auto& path) {
                warn_log_stream << path << ", ";
            });
            warn_log_stream << duplicates.back();
            if (num_duplicate_paths > duplicates.size()) {
                warn_log_stream << " (showing unique duplicates)";
            }
        }
    }
}

std::vector<fs::path> get_read_paths(const OptionMap& options, const bool log = true)
{
    using namespace utils;
    std::vector<fs::path> result {};
    if (is_set("reads", options)) {
        auto resolved_paths = resolve_paths(options.at("reads").as<std::vector<fs::path>>(), options);
        append(std::move(resolved_paths), result);
    }
    if (is_set("reads-file", options)) {
        auto paths_to_read_paths = options.at("reads-file").as<std::vector<fs::path>>();
        for (auto& path_to_read_paths : paths_to_read_paths) {
            path_to_read_paths = resolve_path(path_to_read_paths, options);
            if (!fs::exists(path_to_read_paths)) {
                MissingReadPathFile e {path_to_read_paths};
                e.set_location_specified("the command line option '--reads-file'");
                throw e;
            }
            try {
                auto paths = get_resolved_paths_from_file(path_to_read_paths, options);
                if (log && paths.empty()) {
                    logging::WarningLogger log {};
                    stream(log) << "The read path file you specified " << path_to_read_paths
                                << " in the command line option '--reads-file' is empty";
                }
                append(std::move(paths), result);
            } catch (...) {
                MalfordReadPathFile e {path_to_read_paths};
                e.set_location_specified("the command line option '--reads-file'");
                throw e;
            }
        }
    }
    remove_duplicates(result, "read", log);
    return result;
}

unsigned count_read_paths(const OptionMap& options)
{
    return get_read_paths(options, false).size();
}

ReadManager make_read_manager(const OptionMap& options)
{
    auto read_paths = get_read_paths(options);
    const auto max_open_files = as_unsigned("max-open-read-files", options);
    return ReadManager {std::move(read_paths), max_open_files};
}

bool denovo_candidate_variant_discovery_enabled(const OptionMap& options)
{
    return !options.at("disable-denovo-variant-discovery").as<bool>();
}

bool pileup_candidate_variant_generator_enabled(const OptionMap& options)
{
    return denovo_candidate_variant_discovery_enabled(options)
           && !options.at("disable-pileup-candidate-generator").as<bool>();
}

bool repeat_candidate_variant_generator_enabled(const OptionMap& options)
{
    return denovo_candidate_variant_discovery_enabled(options)
           && !options.at("disable-repeat-candidate-generator").as<bool>();
}

bool assembler_candidate_variant_generator_enabled(const OptionMap& options)
{
    return denovo_candidate_variant_discovery_enabled(options)
        && !options.at("disable-assembly-candidate-generator").as<bool>()
        && !is_fast_mode(options);
}

bool read_preprocessing_disabled(const OptionMap& options)
{
    return options.at("disable-read-preprocessing").as<bool>();
}

AlignedRead::BaseQuality get_max_base_quality(const OptionMap& options)
{
    static constexpr AlignedRead::BaseQuality hard_max {125};
    if (options.count("max-base-quality")) {
        return std::min(options.at("max-base-quality").as<int>(), static_cast<int>(hard_max));
    } else {
        return hard_max;
    }
}

bool soft_clip_masking_enabled(const OptionMap& options)
{
    return options.at("mask-soft-clipped-bases").as<bool>();
}

bool adapter_masking_enabled(const OptionMap& options)
{
    return !options.at("disable-adapter-masking").as<bool>();
}

bool overlap_masking_enabled(const OptionMap& options)
{
    return !options.at("disable-overlap-masking").as<bool>();
}

auto make_read_transformers(const ReferenceGenome& reference, const OptionMap& options)
{
    using namespace octopus::readpipe;
    ReadTransformer prefilter_transformer {}, postfilter_transformer {};
    prefilter_transformer.add(CapitaliseBases {});
    prefilter_transformer.add(CapBaseQualities {get_max_base_quality(options)});
    if (!read_preprocessing_disabled(options)) {
        if (is_set("mask-tails", options)) {
            const auto mask_length = static_cast<MaskTail::Length>(options.at("mask-tails").as<int>());
            prefilter_transformer.add(MaskTail {mask_length});
        }
        if (is_set("mask-low-quality-tails", options)) {
            const auto threshold = static_cast<AlignedRead::BaseQuality>(as_unsigned("mask-low-quality-tails", options));
            prefilter_transformer.add(MaskLowQualityTails {threshold});
        }
        if (soft_clip_masking_enabled(options)) {
            const auto boundary_size = as_unsigned("mask-soft-clipped-boundary-bases", options);
            if (boundary_size > 0) {
                if (is_set("soft-clip-mask-threshold", options)) {
                    const auto threshold = static_cast<AlignedRead::BaseQuality>(as_unsigned("soft-clip-mask-threshold", options));
                    prefilter_transformer.add(MaskLowQualitySoftClippedBoundaryBases {boundary_size, threshold});
                } else if (assembler_candidate_variant_generator_enabled(options)) {
                    prefilter_transformer.add(MaskLowQualitySoftClippedBoundaryBases {boundary_size, 3});
                    prefilter_transformer.add(MaskLowAverageQualitySoftClippedTails {10, 5});
                    prefilter_transformer.add(MaskClippedDuplicatedBases {});
                } else {
                    prefilter_transformer.add(MaskSoftClippedBoundraryBases {boundary_size});
                }
            } else {
                if (is_set("soft-clip-mask-threshold", options)) {
                    const auto threshold = static_cast<AlignedRead::BaseQuality>(as_unsigned("soft-clip-mask-threshold", options));
                    prefilter_transformer.add(MaskLowQualitySoftClippedBases {threshold});
                } else if (assembler_candidate_variant_generator_enabled(options)) {
                    prefilter_transformer.add(MaskLowQualitySoftClippedBases {3});
                    prefilter_transformer.add(MaskLowAverageQualitySoftClippedTails {10, 5});
                    prefilter_transformer.add(MaskClippedDuplicatedBases {});
                } else {
                    prefilter_transformer.add(MaskSoftClipped {});
                }
            }
        }
        if (adapter_masking_enabled(options)) {
            prefilter_transformer.add(MaskAdapters {});
            postfilter_transformer.add(MaskTemplateAdapters {});
        }
        if (overlap_masking_enabled(options)) {
            postfilter_transformer.add(MaskStrandOfDuplicatedBases {});
        }
        if (options.at("mask-inverted-soft-clipping").as<bool>()) {
            prefilter_transformer.add(MaskInvertedSoftClippedReadEnds {reference, 10, 500});
        }
        if (options.at("mask-3prime-shifted-soft-clipped-heads").as<bool>()) {
            prefilter_transformer.add(Mask3PrimeShiftedSoftClippedHeads {reference, 10, 500});
        }
        postfilter_transformer.add(ClearAnnotations {});
        prefilter_transformer.shrink_to_fit();
        postfilter_transformer.shrink_to_fit();
    }
    return std::make_pair(std::move(prefilter_transformer), std::move(postfilter_transformer));
}

bool split_long_reads(const OptionMap& options)
{
    return options.at("split-long-reads").as<bool>();
}

boost::optional<AlignedRead::NucleotideSequence::size_type> max_read_length(const OptionMap& options)
{
    if (split_long_reads(options) && is_set("max-read-length", options)) {
        return as_unsigned("max-read-length", options);
    } else {
        return boost::none;
    }
}

bool use_same_read_profile_for_all_samples(const OptionMap& options)
{
    return options.at("use-same-read-profile-for-all-samples").as<bool>();
}

class BadSamTag : public UserError
{
    std::string name_;
    
    std::string do_where() const override
    {
        return "make_read_filterer";
    }
    std::string do_why() const override
    {
        return "The tag specified in no-reads-with-tag " + name_ + " is invalid";
    }
    std::string do_help() const override
    {
        return "Check the SAM specification (https://samtools.github.io/hts-specs/SAMv1.pdf)";
    }
public:
    BadSamTag(std::string name) : name_ {std::move(name)} {}
};

auto make_read_filterer(const OptionMap& options)
{
    using std::make_unique;
    using namespace octopus::readpipe;
    using ReadFilterer = ReadPipe::ReadFilterer;
    
    ReadFilterer result {};
    
    // these filters are mandatory
    result.add(make_unique<HasValidBaseQualities>());
    result.add(make_unique<HasWellFormedCigar>());
    
    if (read_preprocessing_disabled(options)) {
        return result;
    }
    if (!options.at("consider-unmapped-reads").as<bool>()) {
        result.add(make_unique<IsMapped>());
    }
    
    const auto min_mapping_quality = as_unsigned("min-mapping-quality", options);
    const auto min_base_quality    = as_unsigned("good-base-quality", options);
    const auto min_good_bases      = as_unsigned("min-good-bases", options);
    
    if (min_mapping_quality > 0) {
        result.add(make_unique<IsGoodMappingQuality>(min_mapping_quality));
    }
    if (min_base_quality > 0 && min_good_bases > 0) {
        result.add(make_unique<HasSufficientGoodQualityBases>(min_base_quality, min_good_bases));
    }
    if (min_base_quality > 0 && is_set("min-good-base-fraction", options)) {
        auto min_good_base_fraction = options.at("min-good-base-fraction").as<double>();
        result.add(make_unique<HasSufficientGoodBaseFraction>(min_base_quality, min_good_base_fraction));
    }
    if (is_set("min-read-length", options)) {
        result.add(make_unique<IsLong>(as_unsigned("min-read-length", options)));
    } else {
        result.add(make_unique<IsLong>(0));
    }
    if (is_set("max-read-length", options) && !split_long_reads(options)) {
        result.add(make_unique<IsShort>(as_unsigned("max-read-length", options)));
    }
    if (!options.at("allow-marked-duplicates").as<bool>()) {
        result.add(make_unique<IsNotMarkedDuplicate>());
    }
    if (!options.at("allow-octopus-duplicates").as<bool>()) {
        using RDDP = ReadDeduplicationDetectionPolicy;
        const auto duplicate_detection_policy = options.at("duplicate-read-detection-policy").as<RDDP>();
        switch (duplicate_detection_policy) {
            case RDDP::relaxed: {
                result.add(make_unique<IsNotDuplicate<ReadFilterer::ReadIterator, FivePrimeAndCigarDuplicateDefinition>>());
                break;
            }
            case RDDP::aggressive: {
                result.add(make_unique<IsNotDuplicate<ReadFilterer::ReadIterator, FivePrimeDuplicateDefinition>>());
                break;
            }
        }
    }
    if (!options.at("allow-qc-fails").as<bool>()) {
        result.add(make_unique<IsNotMarkedQcFail>());
    }
    if (is_set("no-reads-with-tag", options)) {
        const auto filter_tags = options.at("no-reads-with-tag").as<std::vector<SamTag>>();
        for (const auto& filter_tag : filter_tags) {
            const auto tag = check_read_tag(filter_tag.tag);
            if (!tag) throw BadSamTag {filter_tag.tag};
            boost::optional<AlignedRead::Annotation> annotation {};
            if (filter_tag.value) annotation = *filter_tag.value;
            result.add(make_unique<NotHasTag>(*tag, std::move(annotation)));
        }
    }
    if (!options.at("allow-secondary-alignments").as<bool>()) {
        result.add(make_unique<IsNotSecondaryAlignment>());
    }
    if (!options.at("allow-supplementary-alignments").as<bool>()) {
        result.add(make_unique<IsNotSupplementaryAlignment>());
    }
    if (options.at("no-reads-with-unmapped-segments").as<bool>()) {
        result.add(make_unique<IsNextSegmentMapped>());
        result.add(make_unique<IsProperTemplate>());
    }
    if (options.at("no-reads-with-distant-segments").as<bool>()) {
        result.add(make_unique<IsLocalTemplate>());
    }
    if (options.at("no-adapter-contaminated-reads").as<bool>()) {
        result.add(make_unique<IsNotContaminated>());
    }
    const auto max_decoy_supplementary_mq = as_unsigned("max-decoy-supplementary-alignment-mapping-quality", options);
    if (max_decoy_supplementary_mq > 0) {
        result.add(make_unique<NoDecoySupplementaryAlignments>(max_decoy_supplementary_mq));
    } else {
        result.add(make_unique<NoDecoySupplementaryAlignments>());
    }
    const auto max_unplaced_supplementary_mq = as_unsigned("max-unplaced-supplementary-alignment-mapping-quality", options);
    if (max_unplaced_supplementary_mq > 0) {
        result.add(make_unique<NoUnplacedSupplementaryAlignments>(max_unplaced_supplementary_mq));
    } else {
        result.add(make_unique<NoUnplacedSupplementaryAlignments>());
    }
    const auto max_unlocalized_supplementary_mq = as_unsigned("max-unlocalized-supplementary-alignment-mapping-quality", options);
    if (max_unlocalized_supplementary_mq > 0) {
        result.add(make_unique<NoUnlocalizedSupplementaryAlignments>(max_unlocalized_supplementary_mq));
    } else {
        result.add(make_unique<NoUnlocalizedSupplementaryAlignments>());
    }
    result.shrink_to_fit();
    return result;
}

bool is_downsampling_enabled(const OptionMap& options)
{
    return !read_preprocessing_disabled(options) && !options.at("disable-downsampling").as<bool>();
}

boost::optional<readpipe::Downsampler> make_downsampler(const OptionMap& options)
{
    if (is_downsampling_enabled(options)) {
        using namespace octopus::readpipe;
        const auto max_coverage    = as_unsigned("downsample-above", options);
        const auto target_coverage = as_unsigned("downsample-target", options);
        return Downsampler {max_coverage, target_coverage};
    }
    return boost::none;
}

ReadPipe make_read_pipe(ReadManager& read_manager, const ReferenceGenome& reference, std::vector<SampleName> samples, const OptionMap& options)
{
    auto transformers = make_read_transformers(reference, options);
    if (transformers.second.num_transforms() > 0) {
        if (split_long_reads(options)) {
            return ReadPipe {read_manager, as_unsigned("max-read-length", options), std::move(transformers.first),
                             make_read_filterer(options),std::move(transformers.second), make_downsampler(options),
                             std::move(samples)};
        } else {
            return ReadPipe {read_manager, std::move(transformers.first), make_read_filterer(options),
                             std::move(transformers.second), make_downsampler(options), std::move(samples)};
        }
    } else {
        return ReadPipe {read_manager, std::move(transformers.first), make_read_filterer(options),
                         make_downsampler(options), std::move(samples)};
    }
}

auto get_default_germline_inclusion_predicate(const OptionMap& options)
{
    return coretools::KnownCopyNumberInclusionPredicate {static_cast<unsigned>(options.at("organism-ploidy").as<int>())};
}

bool is_single_cell_calling(const OptionMap& options)
{
    return options.at("caller").as<std::string>() == "cell";
}

bool is_cancer_calling(const OptionMap& options)
{
    return options.at("caller").as<std::string>() == "cancer"
        || (!is_single_cell_calling(options) && options.count("normal-samples") == 1);
}

bool is_polyclone_calling(const OptionMap& options)
{
    return options.at("caller").as<std::string>() == "polyclone";
}

auto get_min_credible_vaf_probability(const OptionMap& options)
{
    return options.at("min-candidate-credible-vaf-probability").as<float>();
}

auto get_default_somatic_inclusion_predicate(const OptionMap& options, boost::optional<SampleName> normal = boost::none)
{
    const auto min_credible_vaf = options.at("min-credible-somatic-frequency").as<float>();
    const auto min_credible_vaf_probability = get_min_credible_vaf_probability(options);
    if (normal) {
        return coretools::UnknownCopyNumberInclusionPredicate {*normal, min_credible_vaf, min_credible_vaf_probability};
    } else {
        return coretools::UnknownCopyNumberInclusionPredicate {min_credible_vaf, min_credible_vaf_probability};
    }
}

double get_min_clone_vaf(const OptionMap& options)
{
    return options.at("min-clone-frequency").as<float>();
}

auto get_default_polyclone_inclusion_predicate(const OptionMap& options)
{
    const auto min_vaf = get_min_clone_vaf(options) / 2;
    const auto min_vaf_probability = get_min_credible_vaf_probability(options);
    return coretools::UnknownCopyNumberInclusionPredicate {min_vaf, min_vaf_probability};
}

auto get_default_single_cell_inclusion_predicate(const OptionMap& options)
{
    return coretools::CellInclusionPredicate {};
}

coretools::CigarScanner::Options::InclusionPredicate get_candidate_variant_inclusion_predicate(const OptionMap& options) noexcept
{
    if (is_cancer_calling(options)) {
        boost::optional<SampleName> normal{};
        if (is_set("normal-samples", options)) {
            normal = options.at("normal-samples").as<std::vector<SampleName>>().front();
        }
        return get_default_somatic_inclusion_predicate(options, normal);
    } else if (is_polyclone_calling(options)) {
        return get_default_polyclone_inclusion_predicate(options);
    } else if (is_single_cell_calling(options)) {
        return get_default_single_cell_inclusion_predicate(options);
    } else {
        using CVDP = CandidateVariantDiscoveryProtocol;
        if (options.at("variant-discovery-mode").as<CVDP>() == CVDP::illumina) {
            return get_default_germline_inclusion_predicate(options);
        } else {
            return coretools::PacBioInclusionPredicate {};
        }
    }
}

coretools::CigarScanner::Options::MatchPredicate get_candidate_variant_match_predicate(const OptionMap& options)
{
    return coretools::TolerantMatchPredicate {};
}

double get_min_expected_somatic_vaf(const OptionMap& options)
{
    const auto min_credible_frequency = options.at("min-credible-somatic-frequency").as<float>();
    const auto min_expected_frequency = options.at("min-expected-somatic-frequency").as<float>();
    if (std::min(min_credible_frequency, min_expected_frequency) <= 1.0) {
        return std::max(min_credible_frequency, min_expected_frequency);
    } else {
        return std::min(min_credible_frequency, min_expected_frequency);
    }
}

auto get_assembler_region_generator_frequency_trigger(const OptionMap& options)
{
    if (is_cancer_calling(options)) {
        return get_min_expected_somatic_vaf(options);
    } else if (is_polyclone_calling(options)) {
        return get_min_clone_vaf(options);
    } else {
        if (options.at("organism-ploidy").as<int>() < 4) {
            return 0.1;
        } else {
            return 0.05;
        }
    }
}

coretools::LocalReassembler::BubbleScoreSetter
get_assembler_bubble_score_setter(const OptionMap& options) noexcept
{
    using namespace octopus::coretools;
    if (is_cancer_calling(options)) {
        return DepthBasedBubbleScoreSetter {options.at("min-bubble-score").as<double>(),
                                           options.at("min-expected-somatic-frequency").as<float>()};
    } else if (is_polyclone_calling(options)) {
        return DepthBasedBubbleScoreSetter {options.at("min-bubble-score").as<double>(),
                                            options.at("min-clone-frequency").as<float>() / 2};
    } else if (is_single_cell_calling(options)) {
        return DepthBasedBubbleScoreSetter {options.at("min-bubble-score").as<double>(), 0.25};
    } else {
        const auto ploidy = static_cast<unsigned>(options.at("organism-ploidy").as<int>());
        return DepthBasedBubbleScoreSetter {options.at("min-bubble-score").as<double>(), 1. / (10 * ploidy)};
    }
}

boost::optional<double> get_repeat_scanner_min_vaf(const OptionMap& options)
{
    using namespace octopus::coretools;
    if (is_cancer_calling(options)) {
        return options.at("min-credible-somatic-frequency").as<float>() / 2;
    } else if (is_polyclone_calling(options)) {
        return options.at("min-clone-frequency").as<float>() / 2;
    } else if (is_single_cell_calling(options)) {
        return 0.005;
    } else {
        return 0.05;
    }
}

class MissingSourceVariantFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "make_variant_generator_builder";
    }
public:
    MissingSourceVariantFile(fs::path p) : MissingFileError {std::move(p), "source variant"} {};
};

class MissingSourceVariantFileOfPaths : public MissingFileError
{
    std::string do_where() const override
    {
        return "make_variant_generator_builder";
    }
public:
    MissingSourceVariantFileOfPaths(fs::path p) : MissingFileError {std::move(p), "source variant paths"} {};
};

class ConflictingSourceVariantFile : public UserError
{
    std::string do_where() const override
    {
        return "make_variant_generator_builder";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The source variant file you specified " << source_;
        ss << " conflicts with the output file " << output_;
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Specify a unique output file";
    }
    
    fs::path source_, output_;
public:
    ConflictingSourceVariantFile(fs::path source, fs::path output)
    : source_ {std::move(source)}
    , output_ {std::move(output)}
    {}
};

auto get_max_expected_heterozygosity(const OptionMap& options)
{
    const auto snp_heterozygosity = options.at("snp-heterozygosity").as<float>();
    const auto indel_heterozygosity = options.at("indel-heterozygosity").as<float>();
    const auto heterozygosity = snp_heterozygosity + indel_heterozygosity;
    const auto heterozygosity_stdev = options.at("snp-heterozygosity-stdev").as<float>();
    return std::min(static_cast<double>(heterozygosity + 2 * heterozygosity_stdev), 0.9999);
}

auto get_assembler_cycle_tolerance(const OptionMap& options)
{
    using CGT = coretools::LocalReassembler::Options::CyclicGraphTolerance;
    if (options.at("allow-cycles").as<bool>()) {
        using CVDP = CandidateVariantDiscoveryProtocol;
        if (options.at("variant-discovery-mode").as<CVDP>() == CVDP::illumina) {
            return CGT::high;
        } else {
            return CGT::low;
        }
    } else {
        return CGT::none;
    }
}

bool is_big_uncompressed_vcf(const fs::path& vcf_filename)
{
    return vcf_filename.extension().string() == ".vcf" && fs::file_size(vcf_filename) > 500e6;
}

auto make_variant_generator_builder(const OptionMap& options, const boost::optional<const ReadSetProfile&> read_profile)
{
    using namespace coretools;
    
    logging::WarningLogger warning_log {};
    logging::ErrorLogger log {};
    
    VariantGeneratorBuilder result {};
    const bool use_assembler {assembler_candidate_variant_generator_enabled(options)};
    
    if (pileup_candidate_variant_generator_enabled(options)) {
        CigarScanner::Options scanner_options {};
        if (is_set("min-supporting-reads", options)) {
            auto min_support = as_unsigned("min-supporting-reads", options);
            scanner_options.include = coretools::SimpleThresholdInclusionPredicate {min_support};
        } else {
            scanner_options.include = get_candidate_variant_inclusion_predicate(options);
        }
        scanner_options.match = get_candidate_variant_match_predicate(options);
        scanner_options.use_clipped_coverage_tracking = true;
        if (!options.at("force-pileup-candidates").as<bool>()) {
            CigarScanner::Options::MisalignmentParameters misalign_params {};
            misalign_params.max_expected_mutation_rate = get_max_expected_heterozygosity(options);
            misalign_params.snv_threshold = as_unsigned("min-pileup-base-quality", options);
            if (use_assembler) {
                misalign_params.indel_penalty = 1.5;
                misalign_params.clip_penalty = 2;
                misalign_params.min_ln_prob_correctly_aligned = std::log(0.005);
            }
            scanner_options.misalignment_parameters = misalign_params;
        } else {
            scanner_options.misalignment_parameters = boost::none;
        }
        scanner_options.ignore_strand_bias = options.at("allow-strand-biased-candidates").as<bool>();
        result.set_cigar_scanner(std::move(scanner_options));
    }
    if (repeat_candidate_variant_generator_enabled(options)) {
        RepeatScanner::Options repeat_scanner_options {};
        repeat_scanner_options.min_snvs = 1;
        repeat_scanner_options.min_base_quality = 15;
        repeat_scanner_options.min_vaf = get_repeat_scanner_min_vaf(options);
        result.set_repeat_scanner(repeat_scanner_options);
    }
    if (use_assembler) {
        LocalReassembler::Options reassembler_options {};
        const auto kmer_sizes = options.at("kmer-sizes").as<std::vector<int>>();
        reassembler_options.kmer_sizes.assign(std::cbegin(kmer_sizes), std::cend(kmer_sizes));
        if (is_set("assembler-mask-base-quality", options)) {
            reassembler_options.mask_threshold = as_unsigned("assembler-mask-base-quality", options);
        }
        reassembler_options.execution_policy = get_thread_execution_policy(options);
        reassembler_options.num_fallbacks = as_unsigned("max-fallback-kmers", options);
        reassembler_options.fallback_interval_size = as_unsigned("fallback-kmer-gap", options);
        reassembler_options.bin_size = as_unsigned("max-assembly-region-size", options);
        reassembler_options.bin_overlap = as_unsigned("max-assembly-region-overlap", options);
        reassembler_options.min_kmer_observations = as_unsigned("min-kmer-prune", options);
        reassembler_options.max_bubbles = as_unsigned("max-bubbles", options);
        reassembler_options.min_bubble_score = get_assembler_bubble_score_setter(options);
        reassembler_options.max_variant_size = as_unsigned("max-variant-size", options);
        reassembler_options.cycle_tolerance = get_assembler_cycle_tolerance(options);
        reassembler_options.ignore_strand_bias = options.at("allow-strand-biased-candidates").as<bool>();
        result.set_local_reassembler(std::move(reassembler_options));
    }
    if (is_set("source-candidates", options) || is_set("source-candidates-file", options)) {
        const auto output_path = get_output_path(options);
        std::vector<fs::path> source_paths {};
        if (is_set("source-candidates", options)) {
            source_paths = resolve_paths(options.at("source-candidates").as<std::vector<fs::path>>(), options);
        }
        if (is_set("source-candidates-file", options)) {
            auto paths_to_source_paths = options.at("source-candidates-file").as<std::vector<fs::path>>();
            for (auto& path_to_source_paths : paths_to_source_paths) {
                path_to_source_paths = resolve_path(path_to_source_paths, options);
                if (!fs::exists(path_to_source_paths)) {
                    throw MissingSourceVariantFileOfPaths {path_to_source_paths};
                }
                auto file_sources_paths = get_resolved_paths_from_file(path_to_source_paths, options);
                if (file_sources_paths.empty()) {
                    logging::WarningLogger log {};
                    stream(log) << "The source candidate path file you specified " << path_to_source_paths
                                << " in the command line option '--source-candidates-file' is empty";
                }
                utils::append(std::move(file_sources_paths), source_paths);
            }
        }
        remove_duplicates(source_paths, "source variant");
        for (const auto& source_path : source_paths) {
            if (!fs::exists(source_path)) {
                throw MissingSourceVariantFile {source_path};
            }
            if (output_path && source_path == *output_path) {
                throw ConflictingSourceVariantFile {std::move(source_path), *output_path};
            }
            if (is_big_uncompressed_vcf(source_path)) {
                logging::WarningLogger log {};
                stream(log) << "The source candidate VCF you specified " << source_path
                            << " is uncompressed and large - consider compressing and indexing it";
            }
            VcfExtractor::Options vcf_options {};
            vcf_options.max_variant_size = as_unsigned("max-variant-size", options);
            if (is_set("min-source-candidate-quality", options)) {
                vcf_options.min_quality = options.at("min-source-candidate-quality").as<Phred<double>>().score();
            }
            vcf_options.extract_filtered = options.at("use-filtered-source-candidates").as<bool>();
            result.add_vcf_extractor(std::move(source_path), vcf_options);
        }
    }
    if (is_set("regenotype", options)) {
        auto regenotype_path = resolve_path(options.at("regenotype").as<fs::path>(), options);
        if (!fs::exists(regenotype_path)) {
            throw MissingSourceVariantFile {regenotype_path};
        }
        const auto output_path = get_output_path(options);
        if (output_path && regenotype_path == *output_path) {
            throw ConflictingSourceVariantFile {std::move(regenotype_path), *output_path};
        }
        result.add_vcf_extractor(std::move(regenotype_path));
    }
    ActiveRegionGenerator::Options active_region_options {};
    if (is_set("assemble-all", options) && options.at("assemble-all").as<bool>()) {
        active_region_options.assemble_all = true;
    } else {
        AssemblerActiveRegionGenerator::Options assembler_region_options {};
        using TT = AssemblerActiveRegionGenerator::Options::TriggerType;
        assembler_region_options.trigger_types = {TT::indel, TT::structual, TT::clustered};
        assembler_region_options.read_profile = read_profile;
        assembler_region_options.min_expected_mutation_frequency = get_assembler_region_generator_frequency_trigger(options);
        active_region_options.assembler_active_region_generator_options = assembler_region_options;
    }
    result.set_active_region_generator(std::move(active_region_options));
    return result;
}

struct ContigPloidyLess
{
    bool operator()(const ContigPloidy& lhs, const ContigPloidy& rhs) const noexcept
    {
        if (lhs.sample) {
            if (rhs.sample && *lhs.sample != *rhs.sample) {
                return *lhs.sample < *rhs.sample;
            } else {
                return true;
            }
        } else if (rhs.sample) {
            return false;
        }
        return lhs.contig == rhs.contig ? lhs.ploidy < rhs.ploidy : lhs.contig < rhs.contig;
    }
};

struct ContigPloidyEqual
{
    bool operator()(const ContigPloidy& lhs, const ContigPloidy& rhs) const noexcept
    {
        return lhs.sample == rhs.sample && lhs.contig == rhs.contig && lhs.ploidy == rhs.ploidy;
    }
};

struct ContigPloidyAmbiguous
{
    bool operator()(const ContigPloidy& lhs, const ContigPloidy& rhs) const noexcept
    {
        if (lhs.sample && rhs.sample) {
            return *lhs.sample == *rhs.sample && lhs.contig == rhs.contig;
        } else if (!(lhs.sample || rhs.sample)) {
            return lhs.contig == rhs.contig;
        }
        return false;
    }
};

class AmbiguousPloidy : public UserError
{
    std::string do_where() const override
    {
        return "make_caller_factory";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The are contigs with ambiguous ploidy: ";
        for (auto it = std::cbegin(ploidies_), end = std::cend(ploidies_); it != end;) {
            it = std::adjacent_find(it, std::cend(ploidies_), ContigPloidyAmbiguous {});
            if (it != std::cend(ploidies_)) {
                const auto it2 = std::find_if(std::next(it), std::cend(ploidies_),
                                              [=] (const auto& cp) {
                                                  return ContigPloidyAmbiguous{}(*it, cp);
                                              });
                std::ostringstream ss {};
                std::copy(it, it2, std::ostream_iterator<ContigPloidy> {ss, " "});
                it = it2;
            }
        }
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Ensure ploidies are specified only once per sample or per sample contig";
    }
    
    std::vector<ContigPloidy> ploidies_;

public:
    AmbiguousPloidy(std::vector<ContigPloidy> ploidies) : ploidies_ {ploidies} {}
};

void remove_duplicate_ploidies(std::vector<ContigPloidy>& contig_ploidies)
{
    std::sort(std::begin(contig_ploidies), std::end(contig_ploidies), ContigPloidyLess {});
    auto itr = std::unique(std::begin(contig_ploidies), std::end(contig_ploidies), ContigPloidyEqual {});
    contig_ploidies.erase(itr, std::end(contig_ploidies));
}

bool has_ambiguous_ploidies(const std::vector<ContigPloidy>& contig_ploidies)
{
    return std::adjacent_find(std::cbegin(contig_ploidies), std::cend(contig_ploidies),
                              ContigPloidyAmbiguous {}) != std::cend(contig_ploidies);
}

class MissingPloidyFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_ploidy_map";
    }
public:
    MissingPloidyFile(fs::path p) : MissingFileError {std::move(p), "ploidy"} {};
};

PloidyMap get_ploidy_map(const OptionMap& options)
{
    if (options.at("caller").as<std::string>() == "polyclone") {
        return PloidyMap {1};
    }
    std::vector<ContigPloidy> flat_plodies {};
    if (is_set("contig-ploidies-file", options)) {
        const fs::path input_path {options.at("contig-ploidies-file").as<std::string>()};
        const auto resolved_path = resolve_path(input_path, options);
        if (!fs::exists(resolved_path)) {
            throw MissingPloidyFile {input_path};
        }
        std::ifstream file {resolved_path.string()};
        std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                       std::back_inserter(flat_plodies), [] (const Line& line) {
            std::istringstream ss {line.line_data};
            ContigPloidy result {};
            ss >> result;
            return result;
        });
    }
    if (is_set("contig-ploidies", options)) {
        utils::append(options.at("contig-ploidies").as<std::vector<ContigPloidy>>(), flat_plodies);
    }
    remove_duplicate_ploidies(flat_plodies);
    if (has_ambiguous_ploidies(flat_plodies)) {
        throw AmbiguousPloidy {flat_plodies};
    }
    PloidyMap result {as_unsigned("organism-ploidy", options)};
    for (const auto& p : flat_plodies) {
        if (p.sample) {
            result.set(*p.sample, p.contig, p.ploidy);
        } else {
            result.set(p.contig, p.ploidy);
        }
    }
    return result;
}

bool call_sites_only(const OptionMap& options)
{
    return options.at("sites-only").as<bool>();
}

auto get_extension_policy(const OptionMap& options)
{
    using ExtensionPolicy = HaplotypeGenerator::Builder::Policies::Extension;
    switch (options.at("extension-level").as<ExtensionLevel>()) {
        case ExtensionLevel::minimal: return ExtensionPolicy::minimal;
        case ExtensionLevel::conservative: return ExtensionPolicy::conservative;
        case ExtensionLevel::moderate: return ExtensionPolicy::moderate;
        case ExtensionLevel::aggressive: return ExtensionPolicy::aggressive;
        case ExtensionLevel::unlimited: return ExtensionPolicy::unlimited;
        default: return ExtensionPolicy::moderate; // to stop GCC warning
    }
}

auto get_backtrack_policy(const OptionMap& options)
{
    using BacktrackPolicy = HaplotypeGenerator::Builder::Policies::Backtrack;
    switch (options.at("backtrack-level").as<BacktrackLevel>()) {
        case BacktrackLevel::none: return BacktrackPolicy::none;
        case BacktrackLevel::moderate: return BacktrackPolicy::moderate;
        case BacktrackLevel::aggressive: return BacktrackPolicy::aggressive;
        default: return BacktrackPolicy::none; // to stop GCC warning
    }
}

auto get_lagging_policy(const OptionMap& options)
{
    using LaggingPolicy = HaplotypeGenerator::Builder::Policies::Lagging;
    if (is_fast_mode(options)) return LaggingPolicy::none;
    switch (options.at("lagging-level").as<LaggingLevel>()) {
        case LaggingLevel::none: return LaggingPolicy::none;
        case LaggingLevel::conservative: return LaggingPolicy::conservative;
        case LaggingLevel::moderate: return LaggingPolicy::moderate;
        case LaggingLevel::optimistic: return LaggingPolicy::optimistic;
        case LaggingLevel::aggressive: return LaggingPolicy::aggressive;
        default: return LaggingPolicy::none;
    }
}

auto get_max_haplotypes(const OptionMap& options)
{
    if (is_fast_mode(options)) {
        return 50u;
    } else {
        return as_unsigned("max-haplotypes", options);
    }
}

boost::optional<coretools::BadRegionDetector>
make_bad_region_detector(const OptionMap& options, const boost::optional<const ReadSetProfile&>& input_reads_profile)
{
    auto snp_heterozygosity = options.at("snp-heterozygosity").as<float>();
    auto indel_heterozygosity = options.at("indel-heterozygosity").as<float>();
    if (is_cancer_calling(options)) {
        snp_heterozygosity += options.at("somatic-snv-prior").as<float>();
        indel_heterozygosity += options.at("somatic-indel-prior").as<float>();
    }
    const auto heterozygosity = snp_heterozygosity + indel_heterozygosity;
    auto heterozygosity_stdev = options.at("snp-heterozygosity-stdev").as<float>();
    if (is_cancer_calling(options)) {
        heterozygosity_stdev *= (1.0 /  options.at("min-expected-somatic-frequency").as<float>());
    }
    coretools::BadRegionDetector::Parameters params {heterozygosity, heterozygosity_stdev};
    using Tolerance = coretools::BadRegionDetector::Parameters::Tolerance;
    switch (options.at("bad-region-tolerance").as<BadRegionTolerance>()) {
        case BadRegionTolerance::unlimited: return boost::none;
        case BadRegionTolerance::high: params.tolerance = Tolerance::high; break;
        case BadRegionTolerance::normal: params.tolerance = Tolerance::normal; break;
        case BadRegionTolerance::low: params.tolerance = Tolerance::low; break;
    }
    return coretools::BadRegionDetector {params, input_reads_profile};
}

auto get_max_indicator_join_distance() noexcept
{
    return HaplotypeLikelihoodModel{}.pad_requirement();
}

bool allow_flank_scoring(const OptionMap& options)
{
    return !options.at("disable-inactive-flank-scoring").as<bool>() && !is_very_fast_mode(options);
}

auto make_error_model(const OptionMap& options)
{
    const auto& model_label = options.at("sequence-error-model").as<std::string>();
    try {
        return octopus::make_error_model(model_label);
    } catch (const UserError& err) {
        try {
            const auto model_path = resolve_path(model_label, options);
            return octopus::make_error_model(model_path);
        } catch (...) {}
        throw;
    }
}

AlignedRead::MappingQuality calculate_mapping_quality_cap(const OptionMap& options, const boost::optional<const ReadSetProfile&>& read_profile)
{
    constexpr AlignedRead::MappingQuality minimum {60u}; // BWA cap
    if (read_profile) {
        if (read_profile->length_stats.median > 1000) {
            return 240;
        } else if (read_profile->length_stats.median > 200) {
            return 2 * minimum;
        } else {
            return std::max(read_profile->mapping_quality_stats.max, minimum);
        }
    } else {
        return minimum;
    }
}

AlignedRead::MappingQuality calculate_mapping_quality_cap_trigger(const OptionMap& options, const boost::optional<const ReadSetProfile&>& read_profile)
{
    constexpr AlignedRead::MappingQuality minimum {60u}; // BWA cap
    if (read_profile) {
        return std::max(read_profile->mapping_quality_stats.max, minimum);
    } else {
        return minimum;
    }
}

bool model_mapping_quality(const OptionMap& options)
{
    return !options.at("dont-model-mapping-quality").as<bool>();
}

bool use_int_hmm_scores(const OptionMap& options, const boost::optional<const ReadSetProfile&> read_profile)
{
    return options.at("use-wide-hmm-scores").as<bool>();
}

HaplotypeLikelihoodModel
make_calling_haplotype_likelihood_model(const OptionMap& options, const boost::optional<const ReadSetProfile&> read_profile)
{
    auto error_model = make_error_model(options);
    HaplotypeLikelihoodModel::Config config {};
    config.use_mapping_quality = model_mapping_quality(options);
    config.use_flank_state = allow_flank_scoring(options);
    if (config.use_mapping_quality) {
        config.mapping_quality_cap = calculate_mapping_quality_cap(options, read_profile);
        config.mapping_quality_cap_trigger = calculate_mapping_quality_cap_trigger(options, read_profile);
    }
    config.max_indel_error = as_unsigned("max-indel-errors", options);
    config.use_int_scores = use_int_hmm_scores(options, read_profile);
    return HaplotypeLikelihoodModel {std::move(error_model.snv), std::move(error_model.indel), config};
}

bool use_int_hmm_scores_for_realignment(const OptionMap& options, const boost::optional<const ReadSetProfile&> read_profile)
{
    if (options.at("use-wide-hmm-scores").as<bool>()) return true;
    if (options.at("split-long-reads").as<bool>()) return true;
    if (read_profile && read_profile->length_stats.median > 1'000) return true;
    return false;
}

unsigned get_realignment_hmm_max_indel_errors(const OptionMap& options, const boost::optional<const ReadSetProfile&> read_profile)
{
    return as_unsigned("max-indel-errors", options);
}

HaplotypeLikelihoodModel 
make_realignment_haplotype_likelihood_model(const HaplotypeLikelihoodModel& calling_model,
                                            const boost::optional<const ReadSetProfile&> read_profile,
                                            const options::OptionMap& options)
{
    auto result = calling_model;
    auto realignment_config = result.config();
    realignment_config.use_mapping_quality = false;
    realignment_config.use_flank_state = false;
    realignment_config.use_int_scores = use_int_hmm_scores_for_realignment(options, read_profile);
    realignment_config.max_indel_error = get_realignment_hmm_max_indel_errors(options, read_profile);
    result.set(std::move(realignment_config));
    return result;
}

auto get_min_haplotype_flank_pad(const OptionMap& options, const boost::optional<const ReadSetProfile&>& input_reads_profile)
{
    auto model = make_calling_haplotype_likelihood_model(options, input_reads_profile);
    return 2 * model.pad_requirement();
}

auto make_haplotype_generator_builder(const OptionMap& options, const boost::optional<const ReadSetProfile&> input_reads_profile)
{
    const auto lagging_policy    = get_lagging_policy(options);
    const auto max_haplotypes    = get_max_haplotypes(options);
    const auto holdout_limit     = as_unsigned("haplotype-holdout-threshold", options);
    const auto overflow_limit    = as_unsigned("haplotype-overflow", options);
    const auto max_holdout_depth = as_unsigned("max-holdout-depth", options);
    auto result = HaplotypeGenerator::Builder()
        .set_extension_policy(get_extension_policy(options))
        .set_backtrack_policy(get_backtrack_policy(options))
        .set_target_limit(max_haplotypes)
        .set_holdout_limit(holdout_limit)
        .set_overflow_limit(overflow_limit)
        .set_lagging_policy(lagging_policy)
        .set_max_holdout_depth(max_holdout_depth)
        .set_max_indicator_join_distance(get_max_indicator_join_distance())
        .set_min_flank_pad(get_min_haplotype_flank_pad(options, input_reads_profile));
    if (input_reads_profile) {
        result.set_max_allele_distance(1000 * input_reads_profile->length_stats.max);
    }
    return result;
}

boost::optional<Pedigree> read_ped_file(const OptionMap& options)
{
    if (is_set("pedigree", options)) {
        const auto ped_file = resolve_path(options.at("pedigree").as<fs::path>(), options);
        return io::read_pedigree(ped_file);
    } else {
        return boost::none;
    }
}

class BadTrioSampleSet : public UserError
{
    std::string do_where() const override
    {
        return "make_trio";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "Trio calling requires exactly 3 samples but "
           << num_samples_
           << " where provided";
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Ensure only three samples are present; if the read files contain more than"
                " this then explicitly constrain the sample set using the command line option"
                " '--samples'";
    }
    
    std::size_t num_samples_;
    
public:
    BadTrioSampleSet(std::size_t num_samples) : num_samples_ {num_samples} {}
};

class BadTrio : public UserError
{
    std::string do_where() const override
    {
        return "make_trio";
    }
    
    std::string do_why() const override
    {
        return "The given maternal and paternal samples are the same";
    }
    
    std::string do_help() const override
    {
        return "Ensure the sample names given in the command line options"
               " '--maternal-sample' and '--paternal-sample' differ and"
                " refer to valid samples";
    }
};

class BadTrioSamples : public UserError
{
    std::string do_where() const override
    {
        return "make_trio";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        if (mother_ && father_) {
            ss << "Neither of the parent sample names given command line options"
                  " '--maternal-sample' (" << *mother_ << ") and '--paternal-sample' ("
               << *father_ << ") appear in the read sample set";
        } else if (mother_) {
            ss << "The maternal sample name given in the command line option"
                    " '--maternal-sample' (" << *mother_ << ") does not appear in the"
                    " read sample set";
        } else {
            assert(father_);
            ss << "The paternal sample name given in the command line option"
            " '--paternal-sample' (" << *father_  << ") does not appear in the"
            " read sample set";
        }
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Ensure the sample names given in the command line options"
        " '--maternal-sample' and '--paternal-sample' refer to valid samples";
    }
    
    boost::optional<SampleName> mother_, father_;
    
public:
    BadTrioSamples(boost::optional<SampleName> mother, boost::optional<SampleName> father)
    : mother_ {std::move(mother)}
    , father_ {std::move(father)}
    {}
};

auto get_caller_type(const OptionMap& options, const std::vector<SampleName>& samples,
	 				 const boost::optional<Pedigree>& pedigree)
{
    auto result = options.at("caller").as<std::string>();
    if (result == "population" && samples.size() == 1) {
        result = "individual";
    }
    if (is_set("maternal-sample", options) || is_set("paternal-sample", options)
		|| (pedigree && is_trio(samples, *pedigree))) {
        result = "trio";
    }
    if (is_cancer_calling(options)) {
        result = "cancer";
    }
    return result;
}

class BadSampleCount : public UserError
{
    std::string do_where() const override
    {
        return "check_caller";
    }
    
    std::string do_why() const override
    {
        return "The number of samples is not accepted by the chosen caller";
    }
    
    std::string do_help() const override
    {
        return "Check the caller documentation for the required number of samples";
    }
};

template <typename T>
void print_with_quotes(const std::vector<T>& values, std::ostream& os, const char* delim = " ")
{
    if (!values.empty()) {
        std::for_each(std::cbegin(values), std::prev(std::cend(values)), [&] (const auto& value) {
            os << '"' << value << '"' << delim;
        });
        os << '"' << values.back() << '"';
    }
}

class BadNormalSample : public UserError
{
    std::string do_where() const override
    {
        return "check_caller";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The specified normal \"" << given_<< "\""
           << " is not one of the input samples [";
        print_with_quotes(available_, ss);
        ss << "]";
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "Choose one of the input samples as the normal sample";
    }
    SampleName given_;
    std::vector<SampleName> available_;
public:
    BadNormalSample(SampleName given, std::vector<SampleName> available)
    : given_ {std::move(given)}
    , available_ {std::move(available)}
    {}
};

void check_normal_samples(const std::vector<SampleName>& normals, const std::vector<SampleName>& samples)
{
    for (const auto& normal : normals) {
        if (std::find(std::cbegin(samples), std::cend(samples), normal) == std::cend(samples)) {
            throw BadNormalSample {normal, samples};
        }
    }
}

void check_caller(const std::string& caller, const std::vector<SampleName>& samples, const OptionMap& options)
{
    if (caller == "polyclone") {
        if (samples.size() != 1) {
            throw BadSampleCount {};
        }
    } else if (caller == "cancer" || caller == "cell") {
        if (is_set("normal-samples", options)) {
            check_normal_samples(options.at("normal-samples").as<std::vector<std::string>>(), samples);
        }
    }
}

auto get_child_from_trio(std::vector<SampleName> trio, const Pedigree& pedigree)
{
	if (is_parent_of(trio[0], trio[1], pedigree)) return trio[1];
	return is_parent_of(trio[1], trio[0], pedigree) ? trio[0] : trio[2];
}

Trio make_trio(std::vector<SampleName> samples, const Pedigree& pedigree)
{
	return *make_trio(get_child_from_trio(samples, pedigree), pedigree);
}

Trio make_trio(std::vector<SampleName> samples, const OptionMap& options,
			   const boost::optional<Pedigree>& pedigree)
{
	if (pedigree && is_trio(samples, *pedigree)) {
		return make_trio(samples, *pedigree);
	}
    if (samples.size() != 3) {
        throw BadTrioSampleSet {samples.size()};
    }
    auto mother = options.at("maternal-sample").as<SampleName>();
    auto father = options.at("paternal-sample").as<SampleName>();
    if (mother == father) {
        throw BadTrio {};
    }
    std::array<SampleName, 2> parents {mother, father};
    std::vector<SampleName> children {};
    std::sort(std::begin(samples), std::end(samples));
    std::sort(std::begin(parents), std::end(parents));
    assert(std::unique(std::begin(samples), std::end(samples)) == std::end(samples));
    std::set_difference(std::cbegin(samples), std::cend(samples),
                        std::cbegin(parents), std::cend(parents),
                        std::back_inserter(children));
    if (children.size() != 1) {
        boost::optional<SampleName> bad_mother, bad_father;
        if (!std::binary_search(std::cbegin(samples), std::cend(samples), mother)) {
            bad_mother = std::move(mother);
        }
        if (!std::binary_search(std::cbegin(samples), std::cend(samples), father)) {
            bad_father = std::move(father);
        }
        throw BadTrioSamples {std::move(bad_mother), std::move(bad_father)};
    }
    return Trio {
        Trio::Mother {std::move(mother)},
        Trio::Father {std::move(father)},
        Trio::Child  {std::move(children.front())}
    };
}

boost::optional<Pedigree> get_pedigree(const OptionMap& options, const std::vector<SampleName>& samples)
{
    auto result = read_ped_file(options);
    if (!result) {
        if (samples.size() == 3 && is_set("maternal-sample", options) && is_set("paternal-sample", options)) {
            const auto trio = make_trio(samples, options, boost::none);
            result = Pedigree {};
            using Sex = Pedigree::Member::Sex;
            result->add_founder(Pedigree::Member {trio.mother(), Sex::female});
            result->add_founder(Pedigree::Member {trio.father(), Sex::male});
            result->add_descendant(Pedigree::Member {trio.child(), Sex::hermaphroditic}, trio.mother(), trio.father());
        }
    }
    return result;
}

class UnimplementedCaller : public ProgramError
{
    std::string do_where() const override
    {
        return "get_caller_type";
    }
    
    std::string do_why() const override
    {
        return "The " + caller_ + " caller is not yet implemented. Sorry!";
    }
    
    std::string do_help() const override
    {
        return "please wait for updates";
    }
    
    std::string caller_;

public:
    UnimplementedCaller(std::string caller) : caller_ {caller} {}
};

Caller::ModelPosteriorPolicy get_model_posterior_policy(const OptionMap& options)
{
    Caller::ModelPosteriorPolicy result {};
    switch (options.at("model-posterior").as<ModelPosteriorPolicy>()) {
        case ModelPosteriorPolicy::all: result = Caller::ModelPosteriorPolicy::all; break;
        case ModelPosteriorPolicy::off: result = Caller::ModelPosteriorPolicy::off; break;
        case ModelPosteriorPolicy::special: result = Caller::ModelPosteriorPolicy::special; break;
    }
    return result;
}

auto get_normal_contamination_risk(const OptionMap& options)
{
    auto risk = options.at("normal-contamination-risk").as<NormalContaminationRisk>();
    CallerBuilder::NormalContaminationRisk result {};
    switch (risk) {
        case NormalContaminationRisk::high: result = CallerBuilder::NormalContaminationRisk::high; break;
        case NormalContaminationRisk::low: result = CallerBuilder::NormalContaminationRisk::low; break;
    }
    return result;
}

auto get_target_working_memory(const OptionMap& options)
{
    boost::optional<MemoryFootprint> result {};
    if (is_set("target-working-memory", options)) {
        static const MemoryFootprint min_target_memory {*parse_footprint("100M")};
        result = options.at("target-working-memory").as<MemoryFootprint>();
        auto num_threads = get_num_threads(options);
        if (!num_threads) {
            num_threads = std::thread::hardware_concurrency();
        }
        result = MemoryFootprint {std::max(result->bytes() / *num_threads, min_target_memory.bytes())};
    }
    return result;
}

bool is_experimental_caller(const std::string& caller) noexcept
{
    return caller == "population" || caller == "polyclone" || caller == "cell";
}

ReadLinkageType get_read_linkage_type(const OptionMap& options)
{
    switch (options.at("read-linkage").as<ReadLinkage>()) {
        case ReadLinkage::linked: return ReadLinkageType::paired_and_barcode;
        case ReadLinkage::paired: return ReadLinkageType::paired;
        case ReadLinkage::none: return ReadLinkageType::none;
        default: return ReadLinkageType::paired;
    }
}

bool protect_reference_haplotype(const OptionMap& options)
{
    return !options.at("dont-protect-reference-haplotype").as<bool>();
}

boost::optional<std::size_t> get_max_genotypes(const OptionMap& options, const std::string& caller)
{
    if (options.count("max-genotypes") == 1) {
        return as_unsigned("max-genotypes", options);
    } else if (caller == "cancer") {
        return 5'000;
    } else if (caller == "polyclone") {
        return 10'000;
    } else if (as_unsigned("organism-ploidy", options) > 2) {
        return 20'000;
    } else if (is_fast_mode(options)) {
        return 1'000'000;
    } else {
        return boost::none;
    }
}

boost::optional<std::size_t> get_max_genotype_combinations(const OptionMap& options, const std::string& caller)
{
    if (options.count("max-genotype-combinations") == 1) {
        return as_unsigned("max-genotype-combinations", options);
    } else if (caller == "cell") {
        return 1'000;
    } else if (is_fast_mode(options)) {
        return 10'000;
    } else if (caller == "trio") {
        return 100'000;
    } else if (caller == "population") {
        return 100'000;
    } else {
        return boost::none;
    }
}

class TooManyNormalsError : public UserError
{
    std::string do_where() const override
    {
        return "make_caller_factory";
    }
    std::string do_why() const override
    {
        return "Only one normal sample is allowed by the cancer calling model";
    }
    std::string do_help() const override
    {
        return "Only specify one arguement for --normal-sample";
    }
};

auto get_sample_dropout_concentrations(const OptionMap& options)
{
    std::vector<SampleDropoutConcentrationPair> result {};
    if (is_set("sample-dropout-concentrations", options)) {
        result = options.at("sample-dropout-concentrations").as<std::vector<SampleDropoutConcentrationPair>>();
    }
    return result;
}

auto get_phase_detection_policy(const OptionMap& options)
{
    const auto policy = options.at("phasing-policy").as<PhasingPolicy>();
    switch (policy) {
        case PhasingPolicy::conservative: return false;
        case PhasingPolicy::aggressive: return true;
        case PhasingPolicy::automatic:
        // fall through
        default: return as_unsigned("organism-ploidy", options) < 3 && get_read_linkage_type(options) != ReadLinkageType::paired_and_barcode;
    }
}

CallerFactory make_caller_factory(const ReferenceGenome& reference, ReadPipe& read_pipe,
                                  const InputRegionMap& regions, const OptionMap& options,
                                  const boost::optional<const ReadSetProfile&> read_profile)
{
    CallerBuilder vc_builder {reference, read_pipe,
                              make_variant_generator_builder(options, read_profile),
                              make_haplotype_generator_builder(options, read_profile)};
	const auto pedigree = read_ped_file(options);
    const auto caller = get_caller_type(options, read_pipe.samples(), pedigree);
    check_caller(caller, read_pipe.samples(), options);
    vc_builder.set_caller(caller);
    
    if (is_experimental_caller(caller)) {
        logging::WarningLogger log {};
        stream(log) << "The " << caller << " calling model is still in development. Do not use for production work!";
    }
    
    if (options.at("refcall").as<bool>()) {
        emit_in_development_warning("refcall");
        const auto refcall_merge_threshold = options.at("refcall-block-merge-quality").as<Phred<double>>();
        if (refcall_merge_threshold.score() == 0) {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::positional);
        } else {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::blocked);
            vc_builder.set_refcall_merge_block_threshold(refcall_merge_threshold);
        }
        if (is_set("max-refcall-posterior", options)) {
            vc_builder.set_max_refcall_posterior(options.at("max-refcall-posterior").as<Phred<double>>());
        }
    } else {
        vc_builder.set_refcall_type(CallerBuilder::RefCallType::none);
    }
    auto min_variant_posterior = options.at("min-variant-posterior").as<Phred<double>>();
    
    if (is_set("regenotype", options)) {
        if (caller == "cancer") {
            vc_builder.set_min_variant_posterior(min_variant_posterior);
        } else {
            vc_builder.set_min_variant_posterior(Phred<double> {1});
        }
    } else {
        vc_builder.set_min_variant_posterior(min_variant_posterior);
    }
    vc_builder.set_read_linkage(get_read_linkage_type(options));
    vc_builder.set_ploidies(get_ploidy_map(options));
    vc_builder.set_max_haplotypes(get_max_haplotypes(options));
    vc_builder.set_max_genotypes(get_max_genotypes(options, caller));
    vc_builder.set_max_genotype_combinations(get_max_genotype_combinations(options, caller));
    vc_builder.set_haplotype_extension_threshold(options.at("min-protected-haplotype-posterior").as<double>());
    vc_builder.set_reference_haplotype_protection(protect_reference_haplotype(options));
    vc_builder.set_likelihood_model(make_calling_haplotype_likelihood_model(options, read_profile));
    auto min_phase_score = options.at("min-phase-score").as<Phred<double>>();
    vc_builder.set_min_phase_score(min_phase_score);
    vc_builder.set_early_phase_detection_policy(get_phase_detection_policy(options));
    if (!options.at("use-uniform-genotype-priors").as<bool>()) {
        vc_builder.set_snp_heterozygosity(options.at("snp-heterozygosity").as<float>());
        vc_builder.set_indel_heterozygosity(options.at("indel-heterozygosity").as<float>());
    }
    vc_builder.set_model_based_haplotype_dedup(true);
    vc_builder.set_independent_genotype_prior_flag(options.at("use-independent-genotype-priors").as<bool>());
    if (caller == "cancer") {
        if (is_set("normal-samples", options)) {
            const auto& normals = options.at("normal-samples").as<std::vector<std::string>>();
            if (normals.size() > 1) {
                throw TooManyNormalsError {};
            }
            vc_builder.set_normal_sample(normals.front());
        }
        vc_builder.set_max_somatic_haplotypes(as_unsigned("max-somatic-haplotypes", options));
        vc_builder.set_somatic_snv_prior(options.at("somatic-snv-prior").as<float>());
        vc_builder.set_somatic_indel_prior(options.at("somatic-indel-prior").as<float>());
        vc_builder.set_min_expected_somatic_frequency(options.at("min-expected-somatic-frequency").as<float>());
        vc_builder.set_credible_mass(options.at("somatic-credible-mass").as<float>());
        vc_builder.set_min_credible_somatic_frequency(options.at("min-credible-somatic-frequency").as<float>());
        auto min_somatic_posterior = options.at("min-somatic-posterior").as<Phred<double>>();
        vc_builder.set_min_somatic_posterior(min_somatic_posterior);
        vc_builder.set_normal_contamination_risk(get_normal_contamination_risk(options));
        vc_builder.set_tumour_germline_concentration(options.at("tumour-germline-concentration").as<float>());
    } else if (caller == "trio") {
        vc_builder.set_trio(make_trio(read_pipe.samples(), options, pedigree));
        vc_builder.set_snv_denovo_prior(options.at("denovo-snv-prior").as<float>());
        vc_builder.set_indel_denovo_prior(options.at("denovo-indel-prior").as<float>());
        vc_builder.set_min_denovo_posterior(options.at("min-denovo-posterior").as<Phred<double>>());
    } else if (caller == "polyclone") {
        vc_builder.set_max_clones(as_unsigned("max-clones", options));
        const double clone_prior = options.at("clone-prior").as<float>();
        vc_builder.set_clonality_prior([clone_prior] (unsigned clonality) { return maths::geometric_pdf(clonality, 1 - clone_prior); });
        vc_builder.set_clone_concentration(options.at("clone-concentration").as<float>());
    } else if (caller == "cell") {
        if (is_set("normal-samples", options)) {
            for (auto sample : options.at("normal-samples").as<std::vector<std::string>>()) {
                vc_builder.add_normal_sample(std::move(sample));
            }
        }
        vc_builder.set_dropout_concentration(options.at("dropout-concentration").as<float>());
        for (const SampleDropoutConcentrationPair& p : get_sample_dropout_concentrations(options)) {
            vc_builder.set_dropout_concentration(p.sample, p.concentration);
        }
        vc_builder.set_phylogeny_concentration(options.at("phylogeny-concentration").as<float>());
        vc_builder.set_somatic_snv_prior(options.at("somatic-snv-prior").as<float>());
        vc_builder.set_somatic_indel_prior(options.at("somatic-indel-prior").as<float>());
        vc_builder.set_max_clones(as_unsigned("max-clones", options));
        vc_builder.set_max_copy_losses(as_unsigned("max-copy-loss", options));
        vc_builder.set_max_copy_gains(as_unsigned("max-copy-gain", options));
        vc_builder.set_somatic_cnv_prior(options.at("somatic-cnv-prior").as<float>());
    }
    vc_builder.set_model_posterior_policy(get_model_posterior_policy(options));
    if (is_set("max-vb-seeds", options)) vc_builder.set_max_vb_seeds(as_unsigned("max-vb-seeds", options));
    if (call_sites_only(options) && !is_call_filtering_requested(options)) {
        vc_builder.set_sites_only();
    }
    const auto target_working_memory = get_target_working_memory(options);
    if (target_working_memory) vc_builder.set_target_memory_footprint(*target_working_memory);
    vc_builder.set_execution_policy(get_thread_execution_policy(options));
    auto bad_region_detector = make_bad_region_detector(options, read_profile);
    if (bad_region_detector) {
        vc_builder.set_bad_region_detector(std::move(*bad_region_detector));
    }
    return CallerFactory {std::move(vc_builder)};
}

bool is_call_filtering_requested(const OptionMap& options) noexcept
{
    return !options.at("disable-call-filtering").as<bool>() || options.count("annotations") > 0;
}

std::string get_germline_filter_expression(const OptionMap& options)
{
    return options.at("filter-expression").as<std::string>();
}

std::string get_somatic_filter_expression(const OptionMap& options)
{
    return options.at("somatic-filter-expression").as<std::string>();
}

std::string get_denovo_filter_expression(const OptionMap& options)
{
    return options.at("denovo-filter-expression").as<std::string>();
}

std::string get_refcall_filter_expression(const OptionMap& options)
{
    return options.at("refcall-filter-expression").as<std::string>();
}

bool is_filter_training_mode(const OptionMap& options)
{
    return options.at("disable-call-filtering").as<bool>() && options.count("annotations") > 0;
}

bool all_active_measure_annotations_requested(const OptionMap& options)
{
    if (options.count("annotations") == 1) {
        const auto annotations = options.at("annotations").as<std::vector<std::string>>();
        return annotations.size() == 1 && annotations.front() == "active";
    }  else {
        return false;
    }
}

std::set<std::string> get_requested_measure_annotations(const OptionMap& options)
{
    std::set<std::string> result {};
    if (options.count("annotations") == 1) {
        for (const auto& measure : options.at("annotations").as<std::vector<std::string>>()) {
            result.insert(measure);
        }
    }
    return result;
}

bool is_random_forest_filtering(const OptionMap& options)
{
    return is_set("forest-model", options) || is_set("somatic-forest-model", options);
}

auto get_caller_type(const OptionMap& options, const std::vector<SampleName>& samples)
{
    return get_caller_type(options, samples, get_pedigree(options, samples));
}

bool aggregate_annotations(const OptionMap& options)
{
    return options.at("aggregate-annotations").as<bool>();
}

std::unique_ptr<VariantCallFilterFactory>
make_call_filter_factory(const ReferenceGenome& reference, ReadPipe& read_pipe, const OptionMap& options,
                         boost::optional<fs::path> temp_directory)
{
    if (!temp_directory) temp_directory = "/tmp";
    std::unique_ptr<VariantCallFilterFactory> result {};
    if (is_call_filtering_requested(options)) {
        const auto caller = get_caller_type(options, read_pipe.samples());
        if (is_random_forest_filtering(options)) {
            std::vector<RandomForestFilterFactory::Path> forest_files {};
            std::vector<RandomForestFilterFactory::ForestType> forest_types {};
            if (is_set("forest-model", options)) {
                forest_files.push_back(resolve_path(options.at("forest-model").as<fs::path>(), options));
                forest_types.push_back(RandomForestFilterFactory::ForestType::germline);
            }
            RandomForestFilterFactory::Options forest_options {};
            if (is_set("min-forest-quality", options)) forest_options.min_forest_quality = options.at("min-forest-quality").as<Phred<double>>();
            if (caller == "cancer") {
                if (is_set("somatic-forest-model", options)) {
                    forest_files.push_back(resolve_path(options.at("somatic-forest-model").as<fs::path>(), options));
                    forest_types.push_back(RandomForestFilterFactory::ForestType::somatic);
                } else if (options.at("somatics-only").as<bool>()) {
                    forest_types.front() = RandomForestFilterFactory::ForestType::somatic;
                } else {
                    logging::WarningLogger log {};
                    log << "Both germline and somatic forests must be provided for random forest cancer variant filtering";
                }
                forest_options.use_somatic_forest_for_refcalls = !options.at("use-germline-forest-for-somatic-normals").as<bool>();
            }
            result = std::make_unique<RandomForestFilterFactory>(std::move(forest_files), std::move(forest_types), *temp_directory, forest_options);
        } else {
            if (is_filter_training_mode(options)) {
                result = std::make_unique<TrainingFilterFactory>(get_requested_measure_annotations(options));
            } else {
                auto germline_filter_expression = get_germline_filter_expression(options);
                if (caller == "cancer") {
                    if (options.at("somatics-only").as<bool>()) {
                        result = std::make_unique<ThresholdFilterFactory>("", get_somatic_filter_expression(options),
                                                                          "", get_refcall_filter_expression(options));
                    } else {
                        result = std::make_unique<ThresholdFilterFactory>("", germline_filter_expression,
                                                                          "", get_somatic_filter_expression(options),
                                                                          "", get_refcall_filter_expression(options));
                    }
                } else if (caller == "trio") {
                    auto denovo_filter_expression = get_denovo_filter_expression(options);
                    if (options.at("denovos-only").as<bool>()) {
                        result = std::make_unique<ThresholdFilterFactory>("", denovo_filter_expression,
                                                                          "", get_refcall_filter_expression(options),
                                                                          true, ThresholdFilterFactory::Type::denovo);
                    } else {
                        result = std::make_unique<ThresholdFilterFactory>("", germline_filter_expression,
                                                                          "", denovo_filter_expression,
                                                                          "", get_refcall_filter_expression(options),
                                                                          ThresholdFilterFactory::Type::denovo);
                    }
                } else {
                    result = std::make_unique<ThresholdFilterFactory>(germline_filter_expression);
                }
            }
        }
        if (result) {
            VariantCallFilter::OutputOptions output_options {};
            output_options.emit_sites_only = call_sites_only(options);
            if (all_active_measure_annotations_requested(options)) {
                output_options.annotate_all_active_measures = true;
            } else {
                auto annotations = get_requested_measure_annotations(options);
                output_options.annotations.insert(std::begin(annotations), std::end(annotations));
            }
            output_options.aggregate_allele_annotations = aggregate_annotations(options);
            result->set_output_options(std::move(output_options));
        }
    }
    return result;
}

bool use_calling_read_pipe_for_call_filtering(const OptionMap& options) noexcept
{
    return options.at("use-preprocessed-reads-for-filtering").as<bool>();
}

bool keep_unfiltered_calls(const OptionMap& options) noexcept
{
    return options.at("keep-unfiltered-calls").as<bool>();
}

ReadPipe make_default_filter_read_pipe(ReadManager& read_manager, std::vector<SampleName> samples, const OptionMap& options)
{
    using std::make_unique;
    using namespace readpipe;
    ReadTransformer transformer {};
    using ReadFilterer = ReadPipe::ReadFilterer;
    ReadFilterer filterer {};
    filterer.add(make_unique<HasValidBaseQualities>());
    filterer.add(make_unique<HasWellFormedCigar>());
    filterer.add(make_unique<IsMapped>());
    filterer.add(make_unique<IsLong>(1));
    if (!options.at("allow-qc-fails").as<bool>()) {
        filterer.add(make_unique<IsNotMarkedQcFail>());
    }
    if (is_set("no-reads-with-tag", options)) {
        const auto filter_tags = options.at("no-reads-with-tag").as<std::vector<SamTag>>();
        for (const auto& filter_tag : filter_tags) {
            const auto tag = check_read_tag(filter_tag.tag);
            if (!tag) throw BadSamTag {filter_tag.tag};
            boost::optional<AlignedRead::Annotation> annotation {};
            if (filter_tag.value) annotation = *filter_tag.value;
            filterer.add(make_unique<NotHasTag>(*tag, std::move(annotation)));
        }
    }
    return ReadPipe {read_manager, std::move(transformer), std::move(filterer), boost::none, std::move(samples)};
}

ReadPipe make_call_filter_read_pipe(ReadManager& read_manager, const ReferenceGenome& reference, std::vector<SampleName> samples, const OptionMap& options)
{
    if (use_calling_read_pipe_for_call_filtering(options)) {
        return make_read_pipe(read_manager, reference, std::move(samples), options);
    } else {
        return make_default_filter_read_pipe(read_manager, std::move(samples), options);
    }
}

boost::optional<fs::path> get_output_path(const OptionMap& options)
{
    if (is_set("output", options)) {
        return resolve_path(options.at("output").as<fs::path>(), options);
    }
    return boost::none;
}

class UnwritableTempDirectory : public SystemError
{
    std::string do_where() const override { return "create_temp_file_directory"; }
    std::string do_why() const override
    {
        std::ostringstream ss{};
        ss << "Failed to create temporary directory " << directory_;
        if (error_) {
            switch (error_->value()) {
                case boost::system::errc::permission_denied: {
                    ss << ": permission denied";
                    break;
                }
                case boost::system::errc::read_only_file_system: {
                    ss << ": read only file system";
                    break;
                }
                case boost::system::errc::not_enough_memory: {
                    ss << ": not enough memory";
                    break;
                }
                case boost::system::errc::filename_too_long: {
                    ss << ": bad path";
                    break;
                }
                case boost::system::errc::io_error: {
                    ss << ": io error";
                    break;
                }
                default: {
                    ss << ": unexpected error (error code - " << *error_ << ")";
                    break;
                }
            }
        }
        return ss.str();
    }
    std::string do_help() const override
    {
        std::ostringstream ss {};
        if (error_) {
            switch(error_->value()) {
                case boost::system::errc::permission_denied: {
                    ss << "Check user has write permissions to " << directory_.parent_path()
                       << " or select another temp directory location";
                    break;
                }
                case boost::system::errc::read_only_file_system: {
                    ss << "Check user has write permissions to " << directory_.parent_path()
                       << " or select another temp directory location";
                    break;
                }
                case boost::system::errc::not_enough_memory: {
                    ss << "Ensure sufficient disk quota is available";
                    break;
                }
                case boost::system::errc::filename_too_long: {
                    ss << "Specify another temp directory name";
                    break;
                }
                default: {
                    ss << "Send a debug report to " << config::BugReport;
                    break;
                }
            }
        }
        return ss.str();
    }
    
    fs::path directory_;
    boost::optional<boost::system::error_code> error_;

public:
    UnwritableTempDirectory(fs::path directory)
    : directory_ {std::move(directory)}
    , error_ {}
    {}
    UnwritableTempDirectory(fs::path directory, boost::system::error_code error)
    : directory_ {std::move(directory)}
    , error_ {error}
    {}
    
    virtual ~UnwritableTempDirectory() override = default;
};

fs::path create_temp_file_directory(const OptionMap& options)
{
    const auto working_directory = get_working_directory(options);
    auto result = working_directory;
    const fs::path temp_dir_base_name {options.at("temp-directory-prefix").as<fs::path>()};
    result /= temp_dir_base_name;
    constexpr unsigned temp_dir_name_count_limit {10'000};
    unsigned temp_dir_counter {2};
    logging::WarningLogger log {};
    boost::system::error_code error_code {};
    while (!fs::create_directory(result, error_code) && temp_dir_counter <= temp_dir_name_count_limit) {
        if (error_code != boost::system::errc::success) {
            // if create_directory returns false and error code is not set then directory already exists
            // https://stackoverflow.com/a/51804969/2970186
            throw UnwritableTempDirectory {result, error_code};
        }
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
        throw UnwritableTempDirectory {result};
    }
    return result;
}

bool keep_temporary_files(const OptionMap& options)
{
    return options.at("keep-temporary-files").as<bool>();
}

boost::optional<fs::path> filter_request(const OptionMap& options)
{
    if (is_call_filtering_requested(options) && is_set("filter-vcf", options)) {
        return resolve_path(options.at("filter-vcf").as<fs::path>(), options);
    }
    return boost::none;
}

bool annotate_filter_output(const OptionMap& options)
{
    return is_set("annotate-filtered-calls", options);
}

boost::optional<fs::path> bamout_request(const OptionMap& options)
{
    if (is_set("bamout", options)) {
        return resolve_path(options.at("bamout").as<fs::path>(), options);
    }
    return boost::none;
}

bool full_bamouts_requested(const OptionMap& options)
{
    return options.at("bamout-type").as<RealignedBAMType>() == RealignedBAMType::full;
}

unsigned max_open_read_files(const OptionMap& options)
{
    return 2 * std::min(as_unsigned("max-open-read-files", options), count_read_paths(options));
}

unsigned estimate_max_open_files(const OptionMap& options)
{
    unsigned result {0};
    result += max_open_read_files(options);
    if (get_output_path(options)) result += 2;
    result += is_debug_mode(options);
    result += is_trace_mode(options);
    result += is_call_filtering_requested(options);
    return result;
}

boost::optional<fs::path> data_profile_request(const OptionMap& options)
{
    if (is_set("data-profile", options)) {
        return resolve_path(options.at("data-profile").as<fs::path>(), options);
    }
    return boost::none;
}

} // namespace options
} // namespace octopus
