//
//  program_options.hpp
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__program_options__
#define __Octopus__program_options__

#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <memory>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "common.hpp"
#include "genomic_region.hpp"
#include "read_filter.hpp"
#include "read_utils.hpp"
#include "variant_caller.hpp"

class ReferenceGenome;
class ReadManager;
struct Downsampler;
class ReadTransform;
class CandidateVariantGenerator;
class VcfWriter;

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace Octopus
{
    namespace Options
    {
    boost::program_options::variables_map parse_options(int argc, const char** argv);
    
    unsigned get_max_threads(const po::variables_map& options);
    
    size_t get_memory_quota(const po::variables_map& options);
    
    ReferenceGenome get_reference(const po::variables_map& options);
    
    SearchRegions get_search_regions(const po::variables_map& options, const ReferenceGenome& reference);
    
    std::vector<SampleIdType> get_samples(const po::variables_map& options);
    
    std::vector<fs::path> get_read_paths(const po::variables_map& options);
    
    ReadManager get_read_manager(const po::variables_map& options);
    
    ReadFilter<ReadContainer::const_iterator> get_read_filter(const po::variables_map& options);
    
    Downsampler<SampleIdType> get_downsampler(const po::variables_map& options);
        
    ReadTransform get_read_transformer(const po::variables_map& options);
    
    CandidateVariantGenerator get_candidate_generator(const po::variables_map& options, ReferenceGenome& reference);
    
    std::unique_ptr<VariantCaller> get_variant_caller(const po::variables_map& options, ReferenceGenome& reference,
                                                      CandidateVariantGenerator& candidate_generator,
                                                      const GenomicRegion::StringType& contig);
    
    std::unique_ptr<VariantCaller> get_variant_caller(const po::variables_map& options, ReferenceGenome& reference,
                                                      CandidateVariantGenerator& candidate_generator);
    
    VcfWriter get_output_vcf(const po::variables_map& options);
    
    } // namespace Options
} // namespace Octopus

#endif /* defined(__Octopus__program_options__) */
