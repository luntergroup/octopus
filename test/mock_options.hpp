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
#include <boost/optional.hpp>

#include "program_options.hpp"
#include "test_common.hpp"

namespace po = boost::program_options;

inline boost::optional<po::variables_map> get_basic_mock_options()
{
    const char *argv[] = {"octopus",
        "--reference", human_reference_fasta.c_str(),
        //"--reference", ecoli_reference_fasta.c_str(),
        
        //"--reads", HG00101.c_str(),
        "--reads", NA12878_low_coverage.c_str(),
        //"--reads", NA12891_high_coverage.c_str(),
        
        //"--reads", NA12878_high_coverage.c_str(),
        //"--reads", NA12878_low_coverage.c_str(), HG00101.c_str(), HG00102.c_str(), HG00103.c_str(),
        //"--reads", NA12878_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12891_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12878_high_coverage.c_str(), NA12891_high_coverage.c_str(),
        //"--reads", ecoli_bam.string().c_str(),
        
        //"--reads", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.NORMAL.7x.compare.bam", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.7x.n25t65s10.bam",
        
        //"--model", "cancer", // default "population"
        "--normal-sample", "HCC1143.NORMAL.30x.compare", // for cancer model
        
        //"--ploidy", "2",
        "--contig-ploidies", "MT=1", "Y=1",
        
        //"--make-blocked-refcalls",
        //"--make-positional-refcalls",
        
        /* input regions */
        
        // for population
        
        //"--regions", "6:93,705,800-93,706,166", // NA12878_low_coverage no phase
        //"--regions", "6:58,605,652-58,606,061",   // NA12878_low_coverage partial phase
        //"--regions", "6:58,605,687-58,605,779",   // NA12878_low_coverage phase strong
        //"--regions", "3:108,803,741-108,803,854", // NA12878_low_coverage phase weak
        
        //"--regions", "6:22,877,929-22,878,012", // NA12878_low_coverage HMM error
        //"--regions", "6:29,915,924-29,916,412",
        //"--regions", "6:144,712,021-144,712,273",
        //"--regions", "6:89,236,310-89,237,082", // NA12878_high_coverage very nice phasing test
        //"--regions", "6:89,236,734-89,236,784",
        
        "-L", "6:29,690,822-29,717,432", // HLA-F
//        "-L", "6:29,794,350-29,799,548", // HLA-G
//        "-L", "6:29,854,439-29,860,076", // HLA-H
//        "-L", "6:29,909,260-29,914,705", // HLA-A
//        "-L", "6:31,320,151-31,326,581", // HLA-B
//        "-L", "6:31,235,509-31,241,047", // HLA-C
//        "-L", "6:30,456,690-30,462,845", // HLA-E
//        "-L", "6:33,030,330-33,060,358", // HLA-DPA1
//        "-L", "6:33,079,213-33,098,080", // HLA-DPB2
//        "-L", "6:32,971,604-32,978,159", // HLA-DOA
//        "-L", "6:32,472,179-32,786,639", // HLA-rest
        
        //"-L", "11:67,447,414-67,559,102",
        
        //"--regions", "21:11,062,774-11,062,920",
        //"--regions", "21:11062880-11063100",
        //"--regions", "21:11,063,185-11,063,327", // weird haplotpes
        
        //"--regions", "3:47,251,793-47,251,839", // NA12878_high_coverage interesting indels
        //"--regions", "3:47,251,808-47,251,847",
        
        //"-L", "4:40,436,510-40,436,694", // NA12878_low_coverage odd genotype call (should be 1/1)
        
        //"-L", "Y:13447283-13447483", // NA12878_low_coverage too many haplotypes
        
        // for cancer
        //"--regions", "6:52,873,970-52,882,816",
        //"--regions", "5:76,747,066-76,747,106", // not a reversion
        //"--regions", "5:76,781,703-76,781,743", // not a reversion
        //"--regions", "5:76,785,333-76,785,478", // HMM error
        
        //"--skip-regions", "1:1,000,000-2,000,000", "1:1,500,000-10,000,000",
        //"--skip-regions-file", human_skip_regions.c_str(),
        
        // read filters
        "--min-supporting-reads", "2",
        "--min-mapping-quality", "20",
        "--min-snp-base-quality", "15",
        //"--allow-marked-duplicates",
        //"--allow-octopus-duplicates",
        
        "--trim-soft-clipped",
        //"--tail-trim-size", "3",
        //"--trim-adapters",
        
        "--reference-cache-size", "20000",
        "--no-downsampling",
        //"--downsample-above", "500",
        //"--downsample-target", "100",
        
        //"--no-candidates-from-alignments",
        //"--candidates-from-assembler",
        //"--kmer-size", "5",
        //"--candidates-from-source", sample_vcf.c_str(),
        //"--regenotype",
        
        "--min-variant-posterior", "5",
        "--min-refcall-posterior", "1",
        "--min-somatic-posterior", "2",
        
        "--somatics-only",
        
        //"--candidates-from-source", "~/test.bcf",
        //"--regenotype",
        
        //"--max-open-read-files", "1",
        
        //"--output", test_out_vcf.c_str(),
        "--output", "~/Genomics/octopus_test/octopus_variants.vcf",
        
        nullptr};
    
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    return Octopus::Options::parse_options(argc, argv);
}

#endif
