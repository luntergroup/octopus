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

#include "utils/path_utils.hpp"
#include "utils/read_stats.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/string_utils.hpp"
#include "utils/append.hpp"
#include "utils/maths.hpp"
#include "basics/phred.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "basics/ploidy_map.hpp"
#include "readpipe/read_pipe_fwd.hpp"
#include "core/tools/coretools.hpp"
#include "core/callers/caller_builder.hpp"
#include "logging/logging.hpp"
#include "io/region/region_parser.hpp"
#include "io/variant/vcf_reader.hpp"
#include "io/variant/vcf_writer.hpp"
#include "exceptions/user_error.hpp"
#include "exceptions/program_error.hpp"
#include "exceptions/system_error.hpp"
#include "exceptions/missing_file_error.hpp"

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
    if (options.count("working-directory") == 1) {
        auto result = expand_user_path(options.at("working-directory").as<fs::path>());
        if (!fs::exists(result) && !fs::is_directory(result)) {
            throw InvalidWorkingDirectory {result};
        }
        return result;
    }
    return fs::current_path();
}

fs::path resolve_path(const fs::path& path, const OptionMap& options)
{
    return ::octopus::resolve_path(path, get_working_directory(options));
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

} // namespace

std::vector<fs::path> extract_paths_from_file(const fs::path& file_path, const OptionMap& options)
{
    std::ifstream file {file_path.string()};
    assert(file.good());
    std::vector<fs::path> result {};
    std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                   std::back_inserter(result), [] (const auto& line) { return line.line_data; });
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
    const bool result {test.is_open()};
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

MemoryFootprint get_target_read_buffer_size(const OptionMap& options)
{
    return options.at("target-read-buffer-footprint").as<MemoryFootprint>();
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

bool is_fast_mode(const OptionMap& options)
{
    return options.at("fast").as<bool>();
}

ReferenceGenome make_reference(const OptionMap& options)
{
    const fs::path input_path {options.at("reference").as<std::string>()};
    auto resolved_path = resolve_path(input_path, options);
    const auto ref_cache_size = options.at("max-reference-cache-footprint").as<MemoryFootprint>().num_bytes();
    try {
        return octopus::make_reference(std::move(resolved_path),
                                       ref_cache_size,
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
        const auto& input_path = options.at("skip-regions-file").as<fs::path>();
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
        const fs::path input_path {options.at("reads-file").as<fs::path>()};
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
    result.register_transform(CapitaliseBases {});
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

bool is_read_filtering_enabled(const OptionMap& options)
{
    return !options.at("disable-read-filtering").as<bool>();
}

auto make_read_filterer(const OptionMap& options)
{
    using std::make_unique;
    using namespace octopus::readpipe;
    using ReadFilterer = ReadPipe::ReadFilterer;
    
    ReadFilterer result {};
    
    // these filters are mandatory
    result.add(make_unique<HasValidBaseQualities>());
    result.add(make_unique<HasWellFormedCigar>());
    
    if (!is_read_filtering_enabled(options)) {
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
    if (min_base_quality > 0 && options.count("min-good-base-fraction") == 1) {
        auto min_good_base_fraction = options.at("min-good-base-fraction").as<double>();
        result.add(make_unique<HasSufficientGoodBaseFraction>(min_base_quality, min_good_base_fraction));
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

bool is_downsampling_enabled(const OptionMap& options)
{
    return is_read_filtering_enabled(options) && !options.at("disable-downsampling").as<bool>();
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

auto get_default_inclusion_predicate()
{
    static const auto sum = [] (const std::vector<unsigned>& observed_qualities) noexcept {
        return std::accumulate(std::cbegin(observed_qualities), std::cend(observed_qualities), 0);
    };
    static const auto erase_low = [] (std::vector<unsigned>& observed_qualities, const unsigned min) {
        observed_qualities.erase(std::remove_if(std::begin(observed_qualities), std::end(observed_qualities),
                                                [=] (const auto q) { return q < min; }),
                                 std::end(observed_qualities));
    };
    static const auto partial_sort = [] (std::vector<unsigned>& observed_qualities, unsigned n) {
        std::partial_sort(std::begin(observed_qualities),
                          std::next(std::begin(observed_qualities), 2),
                          std::end(observed_qualities),
                          std::greater<> {});
    };
    return [] (const Variant& v, const unsigned depth, std::vector<unsigned>& observed_qualities)
    {
        const auto num_observations = observed_qualities.size();
        if (depth < 4) {
            return num_observations > 1 || sum(observed_qualities) >= 20 || is_deletion(v);
        }
        if (is_snp(v)) {
            const auto base_quality_sum = sum(observed_qualities);
            if (depth <= 60) {
                if (num_observations < 2) return false;
                if (base_quality_sum > 100) return true;
                erase_low(observed_qualities, 5);
                if (observed_qualities.size() < 2) return false;
                if (static_cast<double>(observed_qualities.size()) / depth > 0.2) return true;
                partial_sort(observed_qualities, 2);
                return observed_qualities[0] >= 20 && observed_qualities[1] >= 20;
            } else {
                if (num_observations < 3) return false;
                if (base_quality_sum > 150) return true;
                erase_low(observed_qualities, 10);
                if (observed_qualities.size() < 3) return false;
                if (static_cast<double>(observed_qualities.size()) / depth > 0.2) return true;
                partial_sort(observed_qualities, 3);
                return observed_qualities[0] >= 30 && observed_qualities[1] >= 25 && observed_qualities[2] >= 20;
            }
        } else if (is_insertion(v)) {
            if (num_observations == 1 && alt_sequence_size(v) > 8) return false;
            if (depth <= 15) {
                return num_observations > 1;
            } else if (depth <= 30) {
                if (static_cast<double>(num_observations) / depth > 0.45) return true;
                erase_low(observed_qualities, 20);
                return num_observations > 1;
            } else if (depth <= 60) {
                if (num_observations == 1) return false;
                if (static_cast<double>(num_observations) / depth > 0.4) return true;
                erase_low(observed_qualities, 25);
                if (observed_qualities.size() <= 1) return false;
                if (observed_qualities.size() > 2) return true;
                partial_sort(observed_qualities, 2);
                return static_cast<double>(observed_qualities[0]) / alt_sequence_size(v) > 20;
            } else {
                return (num_observations > 2 && static_cast<double>(num_observations) / depth > 0.1)
                       || static_cast<double>(sum(observed_qualities)) / alt_sequence_size(v) > 20;
            }
        } else {
            return num_observations > 1 && static_cast<double>(num_observations) / depth > 0.05;
        }
    };
}

auto get_default_inclusion_predicate(const OptionMap& options) noexcept
{
    using namespace coretools;
    using InclusionPredicate = DynamicCigarScanner::Options::InclusionPredicate;
    const auto caller = options.at("caller").as<std::string>();
    if (caller == "cancer") {
        // TODO: specialise for this case; we need to be careful about low frequency somatics.
        return InclusionPredicate {get_default_inclusion_predicate()};
    } else {
        return InclusionPredicate {get_default_inclusion_predicate()};
    }
}

auto get_default_match_predicate() noexcept
{
    return [] (const Variant& lhs, const Variant& rhs) noexcept
    {
        if (!are_same_type(lhs, rhs) || is_snp(lhs) || is_mnv(lhs)) {
            return lhs == rhs;
        }
        if (is_insertion(lhs) && alt_sequence_size(lhs) == alt_sequence_size(rhs)) {
            const auto& lhs_alt = alt_sequence(lhs);
            const auto& rhs_alt = alt_sequence(rhs);
            return std::count(std::cbegin(lhs_alt), std::cend(lhs_alt), 'N')
                   == std::count(std::cbegin(rhs_alt), std::cend(rhs_alt), 'N');
        }
        return overlaps(lhs, rhs);
    };
}

bool allow_assembler_generation(const OptionMap& options)
{
    return !(is_fast_mode(options) || options.at("disable-assembly-candidate-generator").as<bool>());
}

auto make_variant_generator_builder(const OptionMap& options)
{
    using namespace coretools;
    
    logging::WarningLogger warning_log {};
    logging::ErrorLogger log {};
    
    VariantGeneratorBuilder result {};
    
    if (!options.at("disable-raw-cigar-candidate-generator").as<bool>()) {
        if (options.count("min-supporting-reads") == 1) {
            CigarScanner::Options scanner_options {};
            scanner_options.min_base_quality = as_unsigned("min-base-quality", options);
            scanner_options.min_support = as_unsigned("min-supporting-reads", options);
            if (scanner_options.min_support == 0) {
                warning_log << "The option --min_supporting_reads was set to 0 - assuming this is a typo and setting to 1";
                ++scanner_options.min_support;
            }
            result.set_cigar_scanner(scanner_options);
        } else {
            DynamicCigarScanner::Options scanner_options {
                get_default_inclusion_predicate(options),
                get_default_match_predicate(),
                true
            };
            result.set_dynamic_cigar_scanner(std::move(scanner_options));
        }
    }
    if (allow_assembler_generation(options)) {
        LocalReassembler::Options reassembler_options {};
        const auto kmer_sizes = options.at("kmer-size").as<std::vector<int>>();
        reassembler_options.kmer_sizes.assign(std::cbegin(kmer_sizes), std::cend(kmer_sizes));
        if (options.count("assembler-mask-base-quality") == 1) {
            reassembler_options.mask_threshold = as_unsigned("assembler-mask-base-quality", options);
        }
        reassembler_options.bin_size = as_unsigned("assembler-bin-size", options);
        result.set_local_reassembler(std::move(reassembler_options));
    }
    if (options.count("generate-candidates-from-source") == 1) {
        const auto input_path = options.at("generate-candidates-from-source").as<fs::path>();
        auto resolved_path = resolve_path(input_path, options);
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << input_path
                         << " given in the input option (--generate-candidates-from-source) does not exist";
        }
        result.add_vcf_extractor(std::move(resolved_path));
    }
    if (options.count("regenotype") == 1) {
        auto regenotype_path = options.at("regenotype").as<fs::path>();
        if (options.count("generate-candidates-from-source") == 1) {
            fs::path input_path {options.at("generate-candidates-from-source").as<std::string>()};
            if (regenotype_path != input_path) {
                warning_log << "Running in regenotype mode but given a different source variant file";
            }
            return result;
        }
        auto resolved_path = resolve_path(regenotype_path, options);
        if (!fs::exists(resolved_path)) {
            stream(log) << "The path " << regenotype_path
                        << " given in the input option (--generate-candidates-from-source) does not exist";
        }
        result.add_vcf_extractor(std::move(resolved_path));
    }
    
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
    std::sort(std::begin(contig_ploidies), std::end(contig_ploidies),
              ContigPloidyLess {});
    const auto itr = std::unique(std::begin(contig_ploidies), std::end(contig_ploidies),
                                 ContigPloidyEqual {});
    contig_ploidies.erase(itr, std::end(contig_ploidies));
}

bool has_ambiguous_ploidies(const std::vector<ContigPloidy>& contig_ploidies)
{
    const auto it2 = std::adjacent_find(std::cbegin(contig_ploidies), std::cend(contig_ploidies),
                                        ContigPloidyAmbiguous {});
    return it2 != std::cend(contig_ploidies);
}

class MissingPloidyFile : public MissingFileError
{
    std::string do_where() const override
    {
        return "get_ploidy_map";
    }
public:
    MissingPloidyFile(fs::path p) : MissingFileError {std::move(p), "read path"} {};
};


PloidyMap get_ploidy_map(const OptionMap& options)
{
    std::vector<ContigPloidy> flat_plodies {};
    if (options.count("contig-ploidies-file") == 1) {
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
    if (options.count("contig-ploidies") == 1) {
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

auto get_lagging_policy(const OptionMap& options)
{
    using LaggingPolicy = HaplotypeGenerator::Builder::Policies::Lagging;
    if (is_fast_mode(options)) return LaggingPolicy::none;
    switch (options.at("phasing-level").as<PhasingLevel>()) {
        case PhasingLevel::aggressive: return LaggingPolicy::aggressive;
        case PhasingLevel::conservative: return LaggingPolicy::conservative;
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

auto make_haplotype_generator_builder(const OptionMap& options)
{
    const auto lagging_policy    = get_lagging_policy(options);
    const auto max_haplotypes    = get_max_haplotypes(options);
    const auto holdout_limit     = as_unsigned("haplotype-holdout-threshold", options);
    const auto overflow_limit    = as_unsigned("haplotype-overflow", options);
    const auto max_holdout_depth = as_unsigned("max-holdout-depth", options);
    return HaplotypeGenerator::Builder()
    .set_target_limit(max_haplotypes).set_holdout_limit(holdout_limit).set_overflow_limit(overflow_limit)
    .set_lagging_policy(lagging_policy).set_max_holdout_depth(max_holdout_depth);
}

auto get_caller_type(const OptionMap& options, const std::vector<SampleName>& samples)
{
    // TODO: could think about getting rid of the 'caller' option and just
    // deduce the caller type directly from the options.
    // Will need to report an error if conflicting caller options are given anyway.
    auto result = options.at("caller").as<std::string>();
    if (result == "population" && samples.size() == 1) {
        result = "individual";
    }
    if (options.count("maternal-sample") == 1 || options.count("paternal-sample") == 1) {
        result = "trio";
    }
    if (options.count("normal-sample") == 1) {
        result = "cancer";
    }
    return result;
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
    BadTrioSamples(boost::optional<SampleName> mother,
                   boost::optional<SampleName> father)
    : mother_ {mother}
    , father_ {father}
    {}
};

Trio make_trio(std::vector<SampleName> samples, const OptionMap& options)
{
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
        const auto iter1 = std::find(std::cbegin(children), std::cend(children), mother);
        const auto iter2 = std::find(std::cbegin(children), std::cend(children), father);
        boost::optional<SampleName> mother, father;
        if (iter1 != std::cend(children)) mother = *iter1;
        if (iter2 != std::cend(children)) father = *iter2;
        throw BadTrioSamples {mother, father};
    }
    return Trio {
        Trio::Mother {std::move(mother)},
        Trio::Father {std::move(father)},
        Trio::Child  {std::move(children.front())}
    };
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

bool allow_flank_scoring(const OptionMap& options)
{
    return !(is_fast_mode(options) || options.at("disable-inactive-flank-scoring").as<bool>());
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
    
    const auto caller = get_caller_type(options, read_pipe.samples());
    if (caller == "population") {
        throw UnimplementedCaller {caller};
    }
    vc_builder.set_caller(caller);
    
    if (options.count("refcall") == 1) {
        const auto refcall_type = options.at("refcall").as<RefCallType>();
        if (refcall_type == RefCallType::positional) {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::positional);
        } else {
            vc_builder.set_refcall_type(CallerBuilder::RefCallType::blocked);
        }
        auto min_refcall_posterior = options.at("min-refcall-posterior").as<Phred<double>>();
        vc_builder.set_min_refcall_posterior(min_refcall_posterior);
    } else {
        vc_builder.set_refcall_type(CallerBuilder::RefCallType::none);
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
    vc_builder.set_ploidies(get_ploidy_map(options));
    vc_builder.set_max_haplotypes(get_max_haplotypes(options));
    vc_builder.set_haplotype_extension_threshold(options.at("haplotype-extension-threshold").as<Phred<double>>());
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
        vc_builder.set_trio(make_trio(read_pipe.samples(), options));
        vc_builder.set_denovo_mutation_rate(options.at("denovo-mutation-rate").as<float>());
    }
    
    vc_builder.set_model_filtering(false); // TODO: turn back on when variant filtering is implemented
//    vc_builder.set_model_filtering(!(options.at("disable-call-filtering").as<bool>()
//                                     || options.at("disable-model-filtering").as<bool>()));
    
    if (call_sites_only(options)) {
        vc_builder.set_sites_only();
    }
    vc_builder.set_flank_scoring(allow_flank_scoring(options));
    
    return CallerFactory {std::move(vc_builder)};
}

boost::optional<fs::path> get_final_output_path(const OptionMap& options)
{
    if (options.count("output") == 1) {
        return resolve_path(options.at("output").as<fs::path>(), options);
    }
    return boost::none;
}

VcfWriter make_output_vcf_writer(const OptionMap& options)
{
    auto output = get_final_output_path(options);
    return output ? VcfWriter {std::move(*output)} : VcfWriter {};
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
