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
        
        "--reads", NA12878_low_coverage.c_str(), HG00101.c_str(), HG00102.c_str(), HG00103.c_str(),
        //"--reads", HG00101.c_str(),
        
        //"--reads", NA12878_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12891_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        
        //"--reads", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.NORMAL.7x.compare.chr17.bam", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.7x.n25t65s10.chr17.bam",
        
        //"--reads", NA12878_high_coverage.c_str(), NA12891_high_coverage.c_str(),
        
        //"--reads", ecoli_bam.string().c_str(),
        
        /* model parameters */
        "--model", "cancer", // default "population"
        "--normal-sample", "HG00102", // for cancer model
        
        "--ploidy", "2",
        "--contig-ploidies", "MT=1", "Y=1",
        
        //"--make-blocked-refcalls",
        //"--make-positional-refcalls",
        
        /* input regions */
        //"--regions", "5:157,031,410-157,031,449",
        //"--regions", "11:67503118-67503253",
        
        //"--regions", "2:104,142,854-104,142,925", // population caller fails here
        //"--regions", "2:104,142,897-104,142,936", // but is ok here
        //"--regions", "2:104142870-104142984", // and here
        
        //"--regions", "2:142376817-142376922",
        //"--regions", "20",
        //"--regions", "8:94,786,581-94,786,835", // empty region
        //"--regions", "5:157,031,410-157,031,449", "8:94,786,581-94,786,835",
        //"--regions", "5:157,031,227-157,031,282",
        //"--regions", "5:157,031,211-157,031,269",
        //"--regions", "5:157,030,955-157,031,902",
        //"--regions", "5:157,031,214-157,031,280",
        //"--regions", "2:104,142,525-104,142,742",
        //"--regions", "2:104,141,138-104,141,448", // really nice example of phasing
        //"--regions", "4:141,265,067-141,265,142", // another nice example of phasing
        //"--regions", "7:58,000,994-58,001,147", // very interesting example
        //"--regions", "7:103,614,916-103,615,058",
        //"--regions", "13:28,265,032-28,265,228", // cool phasing region
        //"--regions", "13:96,095,387-96,095,455",
        //"--regions", "11:81,266,010-81,266,122", // all homo-alt
        //"--regions", "11:81,266,084-81,266,123",
        
        //"--regions", "13:36,803,661-36,803,808", // simple insertion alignment error
        //"--regions", "13:36,803,715-36,803,807",
        
        //"--regions", "16:9,299,984-9,300,090", // complex indels
        //"--regions", "13:96,733,039-96,733,124", // complex indels
        //"--regions", "13:96732801-96733349",
        
        //"--regions", "7:122579662-122579817", // complex indels (unverified)
        
        //"--regions", "16:9,378,560-9,378,687", // very complex indel region
        
        //"--regions", "16:62,646,838-62,647,242", // complex phasable insertion and SNPs (in NA12878-HC)
        //"--regions", "16:62,646,885-62,647,013",
        //"--regions", "6:29,921,451-29,921,623",
        //"--regions", "6:34,410,558-34,410,614",
        
        //"--regions", "Y:22,510,914-22,510,973",
        //"--regions", "Y:17,398,830-17,398,917",
        //"--regions", "MT:11,669-11,768", // very high coverage snp
        
        //"--regions", "5:80,465,625-80,465,862", // crazy indel region
        
        //"--regions", "21:22,137,226-22,137,409", // potential cancer spike in (basic)
        //"--regions", "21:22,137,351-22,137,404",
        //"--regions", "21:22,137,271-22,137,310",
        
        //"--regions", "R00000042:686,055-686,094", // ecoli alignment whim
        
        //"--regions", "2:99,042,722-99,043,403", // region should not be phasable, but is getting phased
        //"--regions", "13:33,749,244-33,749,585", // is this a real SNP?
        //"--regions", "13:47,354,553-47,354,592", // real SNP?
        
        //"--regions", "13:36,802,064-36,805,395", // entire region in NA12878_simulated_cancer_basic..
        //"--regions", "13:36,803,643-36,803,828", // just the region containing the somatic
        //"--regions", "13:36,803,666-36,803,705",
        
        //"--regions", "17:1000000-5000000",
        
        //"--regions", "3:36,629,892-36,630,038", // SNP error here due to alignment error in soft clipping in HG00102 + HG00103
        
        //"--regions", "Y:4,313,467-4,314,172",
        
        //"--regions", "4:63,663,404-63,663,639",
        //"--regions", "4:63,663,761-63,663,996",
        "--regions", "4:91,144,739-91,144,822",
        
        // read filters
        "--min-supporting-reads", "1",
        "--min-mapping-quality", "20",
        "--min-snp-base-quality", "20",
        "--no-marked-duplicates",
        "--no-octopus-duplicates",
        
        // read transforms
        "--trim-soft-clipped",
        "--tail-trim-size", "3",
        "--trim-adapters",
        
        "--reference-cache-size", "20000",
        //"--downsample-above", "500",
        //"--downsample-target", "100",
        
        //"--candidates-from-assembler",
        //"--candidates-from-source", "/Users/danielcooke/Genomics/octopus_test/AllVariants.vcf",
        
        "--min-variant-posterior", "5",
        "--min-refcall-posterior", "1",
        "--min-somatic-posterior", "2",
        
        //"--somatics-only",
        
        //"--candidates-from-source", "~/test.bcf",
        //"--regenotype",
        
        //"--max-open-read-files", "1",
        
        "--output", test_out_vcf.c_str(),
        //"--output", "~/test2.bcf",
        
        nullptr};
    
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    return Octopus::Options::parse_options(argc, argv);
}

#endif
