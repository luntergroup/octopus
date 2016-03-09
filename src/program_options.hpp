//
//  program_options.hpp
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__program_options__
#define __Octopus__program_options__

#include <vector>
#include <cstddef>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "downsampler.hpp"
#include "candidate_generator_builder.hpp"
#include "variant_caller_factory.hpp"

class ReferenceGenome;
class ReadManager;
class ReadTransform;
class VcfWriter;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace Octopus
{
namespace Options
{
    enum class ContigOutputOrder
    {
        LexicographicalAscending, LexicographicalDescending,
        ContigSizeAscending, ContigSizeDescending,
        AsInReferenceIndex, AsInReferenceIndexReversed,
        Unspecified
    };
    
    boost::optional<po::variables_map> parse_options(int argc, const char** argv);
    
    bool is_run_command(const po::variables_map& options);
    
    boost::optional<unsigned> get_num_threads(const po::variables_map& options);
    
    std::size_t get_target_read_buffer_size(const po::variables_map& options);
    
    boost::optional<fs::path> get_log_file_name(const po::variables_map& options);
    
    boost::optional<ReferenceGenome> make_reference(const po::variables_map& options);
    
    SearchRegions get_search_regions(const po::variables_map& options,
                                     const ReferenceGenome& reference);
    
    ContigOutputOrder get_contig_output_order(const po::variables_map& options);
    
    std::vector<SampleIdType> get_samples(const po::variables_map& options);
    
    boost::optional<ReadManager> make_read_manager(const po::variables_map& options);
    
    ReadFilterer make_read_filter(const po::variables_map& options);
    
    boost::optional<Downsampler> make_downsampler(const po::variables_map& options);
    
    ReadTransform make_read_transform(const po::variables_map& options);
    
    CandidateGeneratorBuilder make_candidate_generator_builder(const po::variables_map& options,
                                                               const ReferenceGenome& reference);
    
    VariantCallerFactory make_variant_caller_factory(const ReferenceGenome& reference,
                                                     ReadPipe& read_pipe,
                                                     const CandidateGeneratorBuilder& candidate_generator_builder,
                                                     const SearchRegions& regions,
                                                     const po::variables_map& options);
    
    VcfWriter make_output_vcf_writer(const po::variables_map& options);
    
    boost::optional<fs::path> create_temp_file_directory(const po::variables_map& options);
        
} // namespace Options
} // namespace Octopus

#endif /* defined(__Octopus__program_options__) */
