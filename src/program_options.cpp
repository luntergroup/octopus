//
//  program_options.cpp
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "program_options.hpp"

#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <iterator>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <utility>
#include <thread>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "genomic_region.hpp"
#include "reference_genome.hpp"
#include "aligned_read.hpp"
#include "read_manager.hpp"

#include "read_utils.hpp"
#include "read_filters.hpp"
#include "downsampler.hpp"
#include "read_transform.hpp"
#include "read_transformations.hpp"

#include "haplotype_prior_model.hpp"
#include "variant_caller_builder.hpp"

#include "vcf_reader.hpp"
#include "vcf_writer.hpp"

#include "mappable_algorithms.hpp"
#include "string_utils.hpp"

#include "maths.hpp"

namespace Octopus
{
    namespace Options
    {
    
    void conflicting_options(const po::variables_map& vm, const std::string& opt1, const std::string& opt2)
    {
        if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted()) {
            throw std::logic_error(std::string("conflicting options '") + opt1 + "' and '" + opt2 + "'.");
        }
    }
    
    void option_dependency(const po::variables_map& vm, const std::string& for_what,
                           const std::string& required_option)
    {
        if (vm.count(for_what) && !vm[for_what].defaulted())
            if (vm.count(required_option) == 0 || vm[required_option].defaulted()) {
                throw std::logic_error(std::string("option '") + for_what
                                       + "' requires option '" + required_option + "'.");
            }
    }
    
    boost::optional<po::variables_map> parse_options(int argc, const char** argv)
    {
        try {
            po::positional_options_description p;
            
            p.add("command", -1);
            
            po::options_description general("General options");
            general.add_options()
            ("help,h", "produce help message")
            ("version", "output the version number")
            ("verbosity", po::value<unsigned>()->default_value(0),
             "level of logging. Verbosity 0 switches off logging")
            ;
            
            po::options_description backend("Backend options");
            backend.add_options()
            ("max-threads,t", po::value<unsigned>()->default_value(1),
             "maximum number of threads")
            ("memory", po::value<size_t>()->default_value(8000),
             "target memory usage in MB")
            ("reference-cache-size", po::value<size_t>()->default_value(0),
             "the maximum number of bytes that can be used to cache reference sequence")
            ("compress-reads", po::bool_switch()->default_value(false),
             "compress the reads (slower)")
            ("max-open-read-files", po::value<unsigned>()->default_value(200),
             "the maximum number of read files that can be open at one time")
            ;
            
            po::options_description input("Input/output options");
            input.add_options()
            ("reference,R", po::value<std::string>()->required(), "the reference genome file")
            ("reads,I", po::value<std::vector<std::string>>()->multitoken(),
             "space-seperated list of read file paths")
            ("reads-file", po::value<std::string>(), "list of read file paths, one per line")
            ("regions,L", po::value<std::vector<std::string>>()->multitoken(),
             "space-seperated list of one-indexed variant search regions (chrom:begin-end)")
            ("regions-file", po::value<std::string>(),
             "list of one-indexed variant search regions (chrom:begin-end), one per line")
            ("skip-regions", po::value<std::vector<std::string>>()->multitoken(),
             "space-seperated list of one-indexed regions (chrom:begin-end) to skip")
            ("skip-regions-file", po::value<std::string>(),
             "list of one-indexed regions (chrom:begin-end) to skip, one per line")
            ("samples,S", po::value<std::vector<std::string>>()->multitoken(),
             "space-seperated list of sample names to consider")
            ("samples-file", po::value<std::string>(), "list of sample names to consider, one per line")
            ("output,o", po::value<std::string>()->default_value("octopus_calls.vcf"), "write output to file")
            //("log-file", po::value<std::string>(), "path of the output log file")
            ;
            
            po::options_description filters("Read filter options");
            filters.add_options()
            ("allow-unmapped", po::bool_switch()->default_value(false),
             "turns off marked unmapped read filter")
            ("min-mapping-quality", po::value<unsigned>()->default_value(20),
             "reads with smaller mapping quality are ignored")
            ("good-base-quality", po::value<unsigned>()->default_value(20),
             "base quality threshold used by min-good-bases filter")
            ("min-good-base-fraction", po::value<double>(),
             "base quality threshold used by min-good-bases filter")
            ("min-good-bases", po::value<AlignedRead::SizeType>()->default_value(0),
             "minimum number of bases with quality min-base-quality before read is considered")
            ("allow-qc-fails", po::bool_switch()->default_value(false), "filter reads marked as QC failed")
            ("min-read-length", po::value<AlignedRead::SizeType>(), "filter reads shorter than this")
            ("max-read-length", po::value<AlignedRead::SizeType>(), "filter reads longer than this")
            ("allow-marked-duplicates", po::bool_switch()->default_value(false),
             "allows reads marked as duplicate in alignment record")
            ("allow-octopus-duplicates", po::bool_switch()->default_value(false),
             "allows reads considered duplicates by Octopus")
            ("no-secondary-alignmenets", po::bool_switch()->default_value(false),
             "filters reads marked as secondary alignments")
            ("no-supplementary-alignmenets", po::bool_switch()->default_value(false),
             "filters reads marked as supplementary alignments")
            ("no-unmapped-mates", po::bool_switch()->default_value(false),
             "filters reads with unmapped mates")
            ("downsample-above", po::value<unsigned>()->default_value(10000),
             "downsample reads in regions where coverage is over this")
            ("downsample-target", po::value<unsigned>()->default_value(10000),
             "the target coverage for the downsampler")
            ;
            
            po::options_description transforms("Read transform options");
            transforms.add_options()
            ("trim-soft-clipped", po::bool_switch()->default_value(false),
             "trims soft clipped parts of the read")
            ("tail-trim-size", po::value<AlignedRead::SizeType>()->default_value(0),
             "trims this number of bases off the tail of all reads")
            ("trim-adapters", po::bool_switch()->default_value(false),
             "trims any overlapping regions that pass the fragment size")
            ;
            
            po::options_description candidates("Candidate generation options");
            candidates.add_options()
            ("no-candidates-from-alignments", po::bool_switch()->default_value(false),
             "disables candidate variants from aligned reads")
            ("candidates-from-assembler", po::bool_switch()->default_value(false),
             "generate candidate variants with the assembler")
            ("candidates-from-source", po::value<std::string>(),
             "variant file path containing known variants. These variants will automatically become candidates")
            ("regenotype", po::bool_switch()->default_value(false),
             "disables all generators other than source which must be present")
            ("min-snp-base-quality", po::value<unsigned>()->default_value(20),
             "only base changes with quality above this value are considered for snp generation")
            ("min-supporting-reads", po::value<unsigned>()->default_value(1),
             "minimum number of reads that must support a variant if it is to be considered a candidate")
            ("max-variant-size", po::value<AlignedRead::SizeType>()->default_value(100),
             "maximum candidate varaint size from alignmenet CIGAR")
            ("kmer-size", po::value<unsigned>()->default_value(15),
             "k-mer size to use for assembly")
            ("no-cycles", po::bool_switch()->default_value(false), "dissalow cycles in assembly graph")
            ;
            
            po::options_description model("Model options");
            model.add_options()
            ("model", po::value<std::string>()->default_value("population"), "calling model used")
            ("ploidy", po::value<unsigned>()->default_value(2),
             "organism ploidy, all contigs with unspecified ploidy are assumed this ploidy")
            ("contig-ploidies", po::value<std::vector<std::string>>()->multitoken(),
             "ploidy of individual contigs")
            ("contig-ploidies-file", po::value<std::string>(), "list of contig=ploidy pairs, one per line")
            ("normal-sample", po::value<std::string>(), "normal sample used in cancer model")
            ("maternal-sample", po::value<std::string>(), "maternal sample for trio model")
            ("paternal-sample", po::value<std::string>(), "paternal sample for trio model")
            ("transition-prior", po::value<double>()->default_value(0.003),
             "prior probability of a transition snp from the reference")
            ("transversion-prior", po::value<double>()->default_value(0.003),
             "prior probability of a transversion snp from the reference")
            ("insertion-prior", po::value<double>()->default_value(0.003),
             "prior probability of an insertion into the reference")
            ("deletion-prior", po::value<double>()->default_value(0.003),
             "prior probability of a deletion from the reference")
            ("prior-precision", po::value<double>()->default_value(0.003),
             "precision (inverse variance) of the given variant priors")
            ("max-haplotypes", po::value<unsigned>()->default_value(128),
             "the maximum number of haplotypes the model may consider")
            ;
            
            po::options_description calling("Caller options");
            calling.add_options()
            ("min-variant-posterior", po::value<float>()->default_value(20.0),
             "minimum variant call posterior probability (phred scale)")
            ("min-refcall-posterior", po::value<float>()->default_value(10.0),
             "minimum homozygous reference call posterior probability (phred scale)")
            ("min-somatic-posterior", po::value<float>()->default_value(10.0),
             "minimum somaitc mutation call posterior probability (phred scale)")
            ("make-positional-refcalls", po::bool_switch()->default_value(false),
             "caller will output positional REFCALLs")
            ("make-blocked-refcalls", po::bool_switch()->default_value(false),
             "caller will output blocked REFCALLs")
            ("somatics-only", po::bool_switch()->default_value(false),
             "only output somatic calls (for somatic calling models only)")
            ;
            
            po::options_description all("Allowed options");
            all.add(general).add(backend).add(input).add(filters).add(transforms).add(candidates).add(model).add(calling);
            
            po::variables_map vm;
            po::store(po::command_line_parser(argc, argv).options(all).positional(p).run(), vm);
            
            if (vm.count("help")) {
                std::cout << "Usage: octopus <command> [options]" << std::endl;
                std::cout << all << std::endl;
                return vm;
            }
            
            // boost::option cannot handle option dependencies so we must do our own checks
            
            if (vm.count("reads") == 0 && vm.count("reads-file") == 0) {
                throw boost::program_options::required_option {"--reads | --reads-file"};
            }
            
            if (vm.at("model").as<std::string>() == "cancer" && vm.count("normal-sample") == 0) {
                throw std::logic_error {"option normal-sample is required when model=cancer"};
            }
            
            if (vm.at("model").as<std::string>() == "trio" && (vm.count("maternal-sample") == 0 || vm.count("paternal-sample") == 0)) {
                throw std::logic_error {"option maternal-sample and paternal-sample are required when model=trio"};
            }
            
            conflicting_options(vm, "make-positional-refcalls", "make-blocked-refcalls");
            
            po::notify(vm);
            
            return vm;
        } catch (std::exception& e) {
            std::cout << "Option error: " << e.what() << std::endl;
            return boost::none;
        }
    }
    
    struct Line
    {
        std::string line_data;
        
        operator std::string() const
        {
            return line_data;
        }
    };
    
    std::istream& operator>>(std::istream& str, Line& data)
    {
        std::getline(str, data.line_data);
        return str;
    }
    
    boost::optional<std::string> convert_bed_line_to_region_str(const std::string& bed_line)
    {
        const auto tokens = split(bed_line, '\t');
        
        switch (tokens.size()) {
            case 0:
                return boost::none;
            case 1:
                return std::string {tokens[0]};
            case 2:
                // Assume this represents a half range rather than a point
                return std::string {tokens[0] + ':' + tokens[1] + '-'};
            default:
                return std::string {tokens[0] + ':' + tokens[1] + '-' + tokens[2]};
        }
    }
    
    std::function<boost::optional<GenomicRegion>(const std::string&)>
    make_region_parser(const fs::path& region_path, const ReferenceGenome& reference)
    {
        if (region_path.extension().string() == ".bed") {
            return [&] (const std::string& line) -> boost::optional<GenomicRegion>
            {
                auto region_str = convert_bed_line_to_region_str(line);
                if (region_str) {
                    return parse_region(std::move(*region_str), reference);
                }
                return boost::none;
            };
        } else {
            return [&] (const std::string& line) { return parse_region(line, reference); };
        }
    }
    
    std::vector<GenomicRegion> extract_regions_from_file(const fs::path& file_path,
                                                         const ReferenceGenome& reference)
    {
        if (!fs::exists(file_path)) {
            std::cout << "Input error: file does not exist " << file_path.string() << std::endl;
            return {};
        }
        
        std::ifstream file {file_path.string()};
        
        std::vector<boost::optional<GenomicRegion>> parsed_lines {};
        
        std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                       std::back_inserter(parsed_lines), make_region_parser(file_path, reference));
        
        file.close();
        
        std::vector<GenomicRegion> result {};
        result.reserve(parsed_lines.size());
        
        for (auto&& region : parsed_lines) {
            if (region) {
                result.emplace_back(std::move(*region));
            }
        }
        
        result.shrink_to_fit();
        
        return result;
    }
    
    SearchRegions make_search_regions(const std::vector<GenomicRegion>& regions)
    {
        SearchRegions contig_mapped_regions {};
        
        for (const auto& region : regions) {
            contig_mapped_regions[region.get_contig_name()].insert(region);
        }
        
        SearchRegions result {};
        
        for (const auto& contig_regions : contig_mapped_regions) {
            auto covered_contig_regions = get_covered_regions(std::cbegin(contig_regions.second),
                                                              std::cend(contig_regions.second));
            
            result[contig_regions.first].insert(std::make_move_iterator(std::begin(covered_contig_regions)),
                                                std::make_move_iterator(std::end(covered_contig_regions)));
        }
        
        return result;
    }
    
    SearchRegions extract_search_regions(const ReferenceGenome& reference)
    {
        return make_search_regions(get_all_contig_regions(reference));
    }
    
    SearchRegions extract_search_regions(const std::vector<GenomicRegion>& regions,
                                         std::vector<GenomicRegion>& skip_regions)
    {
        auto input_regions = make_search_regions(regions);
        auto skipped = make_search_regions(skip_regions);
        
        SearchRegions result {};
        result.reserve(input_regions.size());
        
        for (const auto& contig_regions : input_regions) {
            result.emplace(contig_regions.first, splice_all(skipped[contig_regions.first], contig_regions.second));
        }
        
        return result;
    }
    
    SearchRegions extract_search_regions(const ReferenceGenome& reference,
                                         std::vector<GenomicRegion>& skip_regions)
    {
        return extract_search_regions(get_all_contig_regions(reference), skip_regions);
    }
    
    boost::optional<fs::path> get_home_dir()
    {
        static const auto result = fs::path(std::getenv("HOME"));
        
        if (fs::is_directory(result)) {
            return result;
        }
        
        return boost::none;
    }
    
    boost::optional<fs::path> expand_user_path(const fs::path& path)
    {
        if (!path.empty() && path.string().front() == '~') {
            if (path.string().size() > 1 && path.string()[1] == '/') {
                const auto home_dir = get_home_dir();
                
                if (home_dir) {
                    return fs::path {home_dir->string() + path.string().substr(1)};
                }
                
                return boost::none;
            } else {
                return boost::none;
            }
        }
        return path;
    }
    
    boost::optional<std::vector<fs::path>> extract_paths_from_file(const fs::path& file_path)
    {
        const auto expanded_path = expand_user_path(file_path);
        
        std::vector<fs::path> result {};
        
        if (!expanded_path) {
            std::cout << "Octopus: could not expand file path " << file_path << std::endl;
            return boost::none;
        }
        
        if (!fs::exists(*expanded_path)) {
            std::cout << "Octopus: could not find file " << expanded_path << std::endl;
            return boost::none;
        }
        
        std::ifstream file {file_path.string()};
        
        if (!file.good()) {
            std::cout << "Octopus: could not open file " << expanded_path << std::endl;
            return boost::none;
        }
        
        std::transform(std::istream_iterator<Line>(file), std::istream_iterator<Line>(),
                       std::back_inserter(result), [] (const Line& line) { return line.line_data; });
        
        const auto it = std::remove_if(std::begin(result), std::end(result),
                                       [] (const auto& path) { return path.empty(); });
        
        result.erase(it, std::end(result));
        
        return result;
    }
    
    struct ExpandedPaths
    {
        ExpandedPaths() = default;
        std::vector<fs::path> good_paths, bad_paths;
        bool has_good() const noexcept { return !good_paths.empty(); }
        bool has_bad() const noexcept { return !bad_paths.empty(); }
    };
    
    ExpandedPaths expand_user_paths(const std::vector<fs::path>& paths)
    {
        ExpandedPaths result {};
        
        result.good_paths.reserve(paths.size());
        
        for (const auto& path : paths) {
            auto expanded_path = expand_user_path(path);
            
            if (expanded_path) {
                result.good_paths.emplace_back(std::move(*expanded_path));
            } else {
                result.bad_paths.emplace_back(std::move(*expanded_path));
            }
        }
        
        result.good_paths.shrink_to_fit();
        result.bad_paths.shrink_to_fit();
        
        return result;
    }
    
    ExpandedPaths expand_user_paths(const std::vector<std::string>& path_strings)
    {
        std::vector<fs::path> paths {std::cbegin(path_strings), std::cend(path_strings)};
        return expand_user_paths(paths);
    }
        
    unsigned get_max_threads(const po::variables_map& options)
    {
        auto result = options.at("max-threads").as<unsigned>();
        
        if (result == 0) {
            result = std::thread::hardware_concurrency(); // just a hint
            if (result == 0) ++result;
        }
        
        return result;
    }
    
    size_t get_memory_quota(const po::variables_map& options)
    {
        return options.at("memory").as<size_t>();
    }
    
    bool is_run_threaded(const po::variables_map& options)
    {
        const auto num_threads = options.at("max-threads").as<unsigned>();
        return num_threads == 0 || num_threads > 1;
    }
    
    boost::optional<ReferenceGenome> make_reference(const po::variables_map& options)
    {
        const auto& path_str = options.at("reference").as<std::string>();
        
        auto path = expand_user_path(path_str);
        
        if (path) {
            if (fs::exists(*path)) {
                const auto ref_cache_size = options.at("reference-cache-size").as<size_t>();
                return ::make_reference(std::move(*path),
                                        static_cast<ReferenceGenome::SizeType>(ref_cache_size),
                                        is_run_threaded(options));
            } else {
                std::cout << "Octopus: could not find reference path " << *path << std::endl;
            }
        } else {
            std::cout << "Octopus: could not expand reference path " << path_str << std::endl;
        }
        
        return boost::none;
    }
    
    template <typename T, typename S>
    std::vector<T>& append(std::vector<T>& target, const std::vector<S>& source)
    {
        target.insert(std::end(target), std::begin(source), std::end(source));
        return target;
    }
    
    template <typename T>
    std::vector<T>& append(std::vector<T>& target, std::vector<T>&& source)
    {
        target.insert(std::end(target),
                      std::make_move_iterator(std::begin(source)),
                      std::make_move_iterator(std::end(source)));
        source.clear();
        source.shrink_to_fit();
        return target;
    }
    
    std::vector<GenomicRegion> parse_regions(const std::vector<std::string>& unparsed_regions,
                                             const ReferenceGenome& reference)
    {
        std::vector<GenomicRegion> result {};
        result.reserve(unparsed_regions.size());
        
        for (const auto& unparsed_region : unparsed_regions) {
            auto parsed_region = parse_region(unparsed_region, reference);
            if (parsed_region) {
                result.emplace_back(std::move(*parsed_region));
            }
        }
        
        return result;
    }
    
    SearchRegions get_search_regions(const po::variables_map& options, const ReferenceGenome& reference)
    {
        std::vector<GenomicRegion> skip_regions {};
        
        if (options.count("skip-regions") == 1) {
            const auto& regions_strings = options.at("skip-regions").as<std::vector<std::string>>();
            auto parsed_regions = parse_regions(regions_strings, reference);
            append(skip_regions, std::move(parsed_regions));
        }
        
        if (options.count("skip-regions-file") == 1) {
            const auto& skip_path = options.at("skip-regions-file").as<std::string>();
            append(skip_regions, extract_regions_from_file(skip_path, reference));
        }
        
        if (options.count("regions") == 0 && options.count("regions-file") == 0) {
            return extract_search_regions(reference, skip_regions);
        } else {
            std::vector<GenomicRegion> input_regions {};
            
            if (options.count("regions") == 1) {
                const auto& region_strings = options.at("regions").as<std::vector<std::string>>();
                auto parsed_regions = parse_regions(region_strings, reference);
                append(input_regions, std::move(parsed_regions));
            }
            
            if (options.count("regions-file") == 1) {
                const auto& regions_path = options.at("regions-file").as<std::string>();
                append(input_regions, extract_regions_from_file(regions_path, reference));
            }
            
            return extract_search_regions(input_regions, skip_regions);
        }
    }
    
    std::vector<SampleIdType> get_samples(const po::variables_map& options)
    {
        std::vector<SampleIdType> result {};
        
        if (options.count("samples") == 1) {
            auto samples = options.at("samples").as<std::vector<std::string>>();
            result.reserve(samples.size());
            std::copy(std::cbegin(samples), std::cend(samples), std::back_inserter(result));
        }
        
        return result;
    }
    
    void print_bad_paths(const std::vector<fs::path>& bad_paths)
    {
        std::cout << "Octopus: the following paths could not be resolved:" << std::endl;
        for (const auto& path : bad_paths) {
            std::cout << "\t" << path.string() << std::endl;
        }
        std::cout << std::endl;
    }
    
    boost::optional<std::vector<fs::path>> get_read_paths(const po::variables_map& options)
    {
        using std::begin; using std::end;
        
        std::vector<fs::path> result {};
        
        if (options.count("reads") == 1) {
            const auto& read_paths = options.at("reads").as<std::vector<std::string>>();
            
            auto expanded_paths = expand_user_paths(read_paths);
            
            if (!expanded_paths.bad_paths.empty()) {
                print_bad_paths(expanded_paths.bad_paths);
                return boost::none;
            }
            
            append(result, std::move(expanded_paths.good_paths));
        }
        
        if (options.count("reads-file") == 1) {
            const auto& read_file_path = options.at("reads-file").as<std::string>();
            
            auto paths = extract_paths_from_file(read_file_path);
            
            if (!paths) return boost::none;
            
            auto expanded_paths = expand_user_paths(*paths);
            
            if (!expanded_paths.bad_paths.empty()) {
                print_bad_paths(expanded_paths.bad_paths);
                return boost::none;
            }
            
            append(result, std::move(expanded_paths.good_paths));
        }
        
        std::sort(begin(result), end(result));
        
        result.erase(std::unique(begin(result), end(result)), end(result));
        
        return result;
    }
    
    boost::optional<ReadManager> make_read_manager(const po::variables_map& options)
    {
        auto read_paths = get_read_paths(options);
        
        if (read_paths) {
            const auto max_open_files = options.at("max-open-read-files").as<unsigned>();
            return ReadManager {std::move(*read_paths), max_open_files};
        }
        
        return boost::none;
    }
    
    ReadFilterer make_read_filter(const po::variables_map& options)
    {
        using QualityType = AlignedRead::QualityType;
        using SizeType    = AlignedRead::SizeType;
        
        ReadFilterer result {};
        
        if (!options.at("allow-unmapped").as<bool>()) {
            result.register_filter(ReadFilters::is_mapped());
        }
        
        auto min_mapping_quality = options.at("min-mapping-quality").as<unsigned>();
        
        if (min_mapping_quality > 0) {
            result.register_filter(ReadFilters::is_good_mapping_quality(min_mapping_quality));
        }
        
        auto min_base_quality = options.at("good-base-quality").as<unsigned>();
        auto min_good_bases = options.at("min-good-bases").as<unsigned>();
        
        if (min_good_bases > 0) {
            result.register_filter(ReadFilters::has_sufficient_good_quality_bases(min_base_quality, min_good_bases));
        }
        
        if (options.count("min-good-base-fraction") == 1) {
            auto min_good_base_fraction =  options.at("min-good-base-fraction").as<double>();
            result.register_filter(ReadFilters::has_good_base_fraction(min_base_quality, min_good_base_fraction));
        }
        
        if (options.count("min-read-length") == 1) {
            result.register_filter(ReadFilters::is_short(options.at("min-read-length").as<SizeType>()));
        }
        
        if (options.count("max-read-length") == 1) {
            result.register_filter(ReadFilters::is_long(options.at("max-read-length").as<SizeType>()));
        }
        
        if (!options.at("allow-marked-duplicates").as<bool>()) {
            result.register_filter(ReadFilters::is_not_marked_duplicate());
        }
        
        if (!options.at("allow-octopus-duplicates").as<bool>()) {
            result.register_filter(ReadFilters::is_not_duplicate());
        }
        
        if (!options.at("allow-qc-fails").as<bool>()) {
            result.register_filter(ReadFilters::is_not_marked_qc_fail());
        }
        
        if (options.at("no-secondary-alignmenets").as<bool>()) {
            result.register_filter(ReadFilters::is_not_secondary_alignment());
        }
        
        if (options.at("no-supplementary-alignmenets").as<bool>()) {
            result.register_filter(ReadFilters::is_not_supplementary_alignment());
        }
        
        if (options.at("no-unmapped-mates").as<bool>()) {
            result.register_filter(ReadFilters::mate_is_mapped());
        }
        
        return result;
    }
    
    Downsampler make_downsampler(const po::variables_map& options)
    {
        auto max_coverage    = options.at("downsample-above").as<unsigned>();
        auto target_coverage = options.at("downsample-target").as<unsigned>();
        
        return Downsampler(max_coverage, target_coverage);
    }
    
    ReadTransform make_read_transform(const po::variables_map& options)
    {
        using SizeType = AlignedRead::SizeType;
        
        ReadTransform result {};
        
        bool trim_soft_clipped = options.at("trim-soft-clipped").as<bool>();
        
        auto tail_trim_size = options.at("tail-trim-size").as<SizeType>();
        
        if (trim_soft_clipped && tail_trim_size > 0) {
            result.register_transform(ReadTransforms::trim_soft_clipped_tails(tail_trim_size));
        } else if (tail_trim_size > 0) {
            result.register_transform(ReadTransforms::trim_tail(tail_trim_size));
        } else if (trim_soft_clipped) {
            result.register_transform(ReadTransforms::trim_soft_clipped());
        }
        
        if (options.at("trim-adapters").as<bool>()) {
            result.register_transform(ReadTransforms::trim_adapters());
        }
        
        return result;
    }
    
    CandidateGeneratorBuilder make_candidate_generator_builder(const po::variables_map& options,
                                                               const ReferenceGenome& reference)
    {
        CandidateGeneratorBuilder result {};
        
        if (options.count("candidates-from-source") == 1) {
            result.add_generator(CandidateGeneratorBuilder::Generator::External);
            
            auto source_path = expand_user_path(options.at("candidates-from-source").as<std::string>());
            
            if (source_path) {
                result.set_variant_source(std::move(*source_path));
            } else {
                std::cout << "Input warning: could not expand candidate source file." << std::endl;
            }
        }
        
        if (options.at("regenotype").as<bool>()) {
            if (options.count("candidates-from-source") == 0) {
                std::cout << "Input warning: source variant file(s) must be present in regenotype mode" << std::endl;
            }
            return result;
        }
        
        if (!options.at("no-candidates-from-alignments").as<bool>()) {
            result.add_generator(CandidateGeneratorBuilder::Generator::Alignment);
            
            result.set_reference(reference);
            result.set_min_snp_base_quality(options.at("min-snp-base-quality").as<unsigned>());
            
            auto min_supporting_reads = options.at("min-supporting-reads").as<unsigned>();
            
            if (min_supporting_reads == 0) ++min_supporting_reads; // probably input error; 0 is meaningless
            
            result.set_min_supporting_reads(min_supporting_reads);
        }
        
        if (options.at("candidates-from-assembler").as<bool>()) {
            result.add_generator(CandidateGeneratorBuilder::Generator::Assembler);
            result.set_kmer_size(options.at("kmer-size").as<unsigned>());
            //auto allow_cycles = !options.at("no-cycles").as<bool>();
        }
        
        result.set_max_variant_size(options.at("max-variant-size").as<CandidateGeneratorBuilder::SizeType>());
        
        return result;
    }
    
//        std::unordered_map<std::string, unsigned> parse_contig_ploidies(const po::variables_map& options)
//        {
//            std::unordered_map<std::string, unsigned> result {};
//            
//            if (options.count("contig-ploidies") == 1) {
//                auto contig_ploidies = options.at("contig-ploidies").as<std::vector<std::string>>();
//                
//                for (const auto& contig_ploidy : contig_ploidies) {
//                    
//                    if (contig_ploidy.find(contig) == 0) {
//                        if (contig_ploidy[contig.size()] != '=') {
//                            throw std::runtime_error {"Could not pass contig-plodies option"};
//                        }
//                        ploidy = static_cast<unsigned>(std::stoul(contig_ploidy.substr(contig.size() + 1)));
//                    }
//                }
//            }
//            
//            if (options.count("contig-ploidies-file") == 1) {
//                
//            }
//            
//            return result;
//        }
        
        unsigned extract_contig_ploidy(const GenomicRegion::ContigNameType& contig,
                                       const po::variables_map& options)
        {
            unsigned result {options.at("ploidy").as<unsigned>()};
            
            if (options.count("contig-ploidies") == 1) {
                auto contig_ploidies = options.at("contig-ploidies").as<std::vector<std::string>>();
                
                for (const auto& contig_ploidy : contig_ploidies) {
                    if (contig_ploidy.find(contig) == 0) {
                        if (contig_ploidy[contig.size()] != '=') {
                            throw std::runtime_error {"Could not pass contig-plodies option"};
                        }
                        result = static_cast<unsigned>(std::stoul(contig_ploidy.substr(contig.size() + 1)));
                    }
                }
            } else if (options.count("contig-ploidies-file") == 1) {
                // TODO: fetch from file
            }
            
            return result;
        }
    
    std::unique_ptr<VariantCaller> make_variant_caller(const po::variables_map& options,
                                                       const ReferenceGenome& reference,
                                                       const CandidateGeneratorBuilder& candidate_generator_builder,
                                                       const GenomicRegion::ContigNameType& contig)
    {
        using Maths::phred_to_probability;
        
        VariantCallerBuilder vc_builder {reference, candidate_generator_builder};
        
        const auto& model = options.at("model").as<std::string>();
        
        vc_builder.set_model(model);
        
        if (options.at("make-positional-refcalls").as<bool>()) {
            vc_builder.set_refcall_type(VariantCaller::RefCallType::Positional);
        } else if (options.at("make-blocked-refcalls").as<bool>()) {
            vc_builder.set_refcall_type(VariantCaller::RefCallType::Blocked);
        }
        
        vc_builder.set_ploidy(extract_contig_ploidy(contig, options));
        
        const auto min_variant_posterior_phred = options.at("min-variant-posterior").as<float>();
        vc_builder.set_min_variant_posterior(phred_to_probability(min_variant_posterior_phred));
        
        const auto min_refcall_posterior_phred = options.at("min-refcall-posterior").as<float>();
        vc_builder.set_min_refcall_posterior(phred_to_probability(min_refcall_posterior_phred));
        
        if (model == "cancer") {
            vc_builder.set_normal_sample(options.at("normal-sample").as<std::string>());
            
            const auto min_somatic_posterior_phred = options.at("min-somatic-posterior").as<float>();
            vc_builder.set_min_somatic_posterior(phred_to_probability(min_somatic_posterior_phred));
            
            if (options.at("somatics-only").as<bool>()) {
                vc_builder.set_somatic_only_calls();
            } else {
                vc_builder.set_somatic_and_variant_calls();
            }
        } else if (model == "trio") {
            vc_builder.set_maternal_sample(options.at("maternal-sample").as<std::string>());
            vc_builder.set_paternal_sample(options.at("paternal-sample").as<std::string>());
        }
        
        return vc_builder.build();
    }
    
    VcfWriter make_output_vcf_writer(const po::variables_map& options)
    {
        auto out_path = expand_user_path(options.at("output").as<std::string>());
        
        if (out_path) {
            return VcfWriter {std::move(*out_path)};
        }
        
        return VcfWriter {};
    }
    
    } // namespace Options
} // namespace Octopus
