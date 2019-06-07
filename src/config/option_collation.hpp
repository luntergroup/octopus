// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef option_collation_hpp
#define option_collation_hpp

#include <vector>
#include <cstddef>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "common.hpp"
#include "option_parser.hpp"
#include "basics/ploidy_map.hpp"
#include "core/callers/caller_factory.hpp"
#include "core/csr/filters/variant_call_filter_factory.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/read/read_manager.hpp"
#include "io/variant/vcf_writer.hpp"
#include "readpipe/read_pipe.hpp"
#include "utils/input_reads_profiler.hpp"
#include "utils/memory_footprint.hpp"

namespace fs = boost::filesystem;

namespace octopus { namespace options {

bool is_run_command(const OptionMap& options);

bool is_debug_mode(const OptionMap& options);
bool is_trace_mode(const OptionMap& options);

boost::optional<fs::path> get_debug_log_file_name(const OptionMap& options);
boost::optional<fs::path> get_trace_log_file_name(const OptionMap& options);

boost::optional<unsigned> get_num_threads(const OptionMap& options);

MemoryFootprint get_target_read_buffer_size(const OptionMap& options);

ReferenceGenome make_reference(const OptionMap& options);

InputRegionMap get_search_regions(const OptionMap& options, const ReferenceGenome& reference);

ContigOutputOrder get_contig_output_order(const OptionMap& options);

bool ignore_unmapped_contigs(const OptionMap& options);

boost::optional<std::vector<SampleName>> get_user_samples(const OptionMap& options);

ReadManager make_read_manager(const OptionMap& options);

ReadPipe make_read_pipe(ReadManager& read_manager, const ReferenceGenome& reference, std::vector<SampleName> samples, const OptionMap& options);

bool call_sites_only(const OptionMap& options);

PloidyMap get_ploidy_map(const OptionMap& options);

boost::optional<Pedigree> get_pedigree(const OptionMap& options, const std::vector<SampleName>& samples);

HaplotypeLikelihoodModel make_haplotype_likelihood_model(const OptionMap& options, const boost::optional<ReadSetProfile>& read_profile);

CallerFactory make_caller_factory(const ReferenceGenome& reference, ReadPipe& read_pipe,
                                  const InputRegionMap& regions, const OptionMap& options,
                                  boost::optional<ReadSetProfile> input_reads_profile = boost::none);

bool is_call_filtering_requested(const OptionMap& options) noexcept;

std::unique_ptr<VariantCallFilterFactory>
make_call_filter_factory(const ReferenceGenome& reference, ReadPipe& read_pipe, const OptionMap& options,
                         boost::optional<fs::path> temp_directory = boost::none);

bool use_calling_read_pipe_for_call_filtering(const OptionMap& options) noexcept;

bool keep_unfiltered_calls(const OptionMap& options) noexcept;

ReadPipe make_call_filter_read_pipe(ReadManager& read_manager, const ReferenceGenome& reference, std::vector<SampleName> samples, const OptionMap& options);

boost::optional<fs::path> get_output_path(const OptionMap& options);

fs::path create_temp_file_directory(const OptionMap& options);

bool is_legacy_vcf_requested(const OptionMap& options);

bool is_filter_training_mode(const OptionMap& options);

boost::optional<fs::path> filter_request(const OptionMap& options);
bool annotate_filter_output(const OptionMap& options);

boost::optional<fs::path> bamout_request(const OptionMap& options);
bool full_bamouts_requested(const OptionMap& options);

unsigned estimate_max_open_files(const OptionMap& options);

boost::optional<fs::path> data_profile_request(const OptionMap& options);

} // namespace options
} // namespace octopus

#endif
