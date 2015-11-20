//
//  mock_options.h
//  Octopus
//
//  Created by Daniel Cooke on 14/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mock_options_h
#define Octopus_mock_options_h

#include <boost/program_options.hpp>

#include "program_options.hpp"
#include "test_common.hpp"

namespace po = boost::program_options;

inline po::variables_map get_basic_mock_options()
{
    const char *argv[] = {"octopus",
        "--reference", human_reference_fasta.c_str(),
        //"--reference", ecoli_reference_fasta.c_str(),
        
        "--reads", NA12878_high_coverage.c_str(),
        //"--reads", NA12878_low_coverage.c_str(),
        
        //"--reads", NA12878_low_coverage.c_str(), HG00101.c_str(), HG00102.c_str(), HG00103.c_str(),
        //"--reads", HG00101.c_str(),
        //"--reads", NA12878_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12891_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.NORMAL.7x.compare.chr17.bam", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.7x.n25t65s10.chr17.bam",
        //"--reads", NA12878_high_coverage.c_str(), NA12891_high_coverage.c_str(),
        //"--reads", ecoli_bam.string().c_str(),
        
        //"--model", "cancer", // default "population"
        "--normal-sample", "HG00102", // for cancer model
        
        //"--ploidy", "2",
        "--contig-ploidies", "MT=1", "Y=1",
        
        //"--make-blocked-refcalls",
        //"--make-positional-refcalls",
        
        //"--regions", "6:29,915,924-29,916,412",
        //"--regions", "6:144,712,021-144,712,273",
        "--regions", "6:89,236,560-89,236,812", // NA12878_high_coverage very nice phasing test
        
        "--min-supporting-reads", "2",
        "--min-mapping-quality", "10",
        "--min-snp-base-quality", "20",
        "--no-marked-duplicates",
        "--no-octopus-duplicates",
        
        "--trim-soft-clipped",
        //"--tail-trim-size", "3",
        //"--trim-adapters",
        
        "--reference-cache-size", "20000",
        //"--downsample-above", "500",
        //"--downsample-target", "100",
        
        //"--no-candidates-from-alignments",
        //"--candidates-from-assembler",
        //"--kmer-size", "5",
        //"--candidates-from-source", "/Users/danielcooke/Genomics/octopus_test/AllVariants.vcf",
        
        "--min-variant-posterior", "5",
        "--min-refcall-posterior", "1",
        "--min-somatic-posterior", "2",
        
        //"--somatics-only",
        
        //"--candidates-from-source", "~/test.bcf",
        //"--regenotype",
        
        //"--max-open-read-files", "1",
        
        //"--output", test_out_vcf.c_str(),
        "--output", "~/test.vcf",
        
        nullptr};
    
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    return Octopus::Options::parse_options(argc, argv);
}

#endif
