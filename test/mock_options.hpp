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
    const char *argv[] = {
        "octopus",
        
        //"--help",
        //"--version",
        
        "--working-directory", "~/",
        
        "--target-read-buffer-size", "0.5",
        //"--threaded",
        
        //"--contig-output-order", "as-in-reference-reversed",
        
        "--reference", human_reference_fasta.c_str(),
        //"--reference", ecoli_reference_fasta.c_str(),
        
        //"--reads", NA12878_low_coverage.c_str(),
        //"--reads", NA12878_high_coverage.c_str(),
        "--reads", "~/Genomics/Illumina/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.chr22.bam",
        
        //"--reads", NA12891_high_coverage.c_str(),
        
        //"--reads", HG00101.c_str(),
        //"--reads", NA12878_low_coverage.c_str(), HG00101.c_str(), HG00102.c_str(), HG00103.c_str(),
        //"--reads", NA12878_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12891_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12878_high_coverage.c_str(), NA12891_high_coverage.c_str(),
        //"--reads", ecoli_bam.string().c_str(),
        
        //"--reads", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.NORMAL.7x.compare.bam", "/Users/dcooke/Genomics/cancer/TCGA/benchmark/HCC1143.7x.n25t65s10.bam",
        
        //"--reads", NA12878_low_coverage_cram.c_str(),
        
        //"--model", "cancer", // default "population"
        "--normal-sample", "NA12878", // for cancer model
        
        //"--organism-ploidy", "3",
        "--contig-ploidies", "MT=1", "Y=1",// "MT=2",
        //"--contig-ploidies-file", contig_ploidies_txt_file.c_str(),
        
        //"--make-blocked-refcalls",
        //"--make-positional-refcalls",
        
        /* input regions */
        
        //"--use-one-based-indexing",
        
        // for population
        
        // NA12878_high_coverage possible errors
        "--regions", "22:37,557,217-37,557,521", // false positive deletion
//        "--regions", "22:37,616,864-37,617,015", // false negative insertion?
//        "--regions", "22:37,622,869-37,622,909", // SNPs or indels?
//        "--regions", "22:37,622,803-37,623,840", // very repetitive region
//        "--regions", "22:37,777,920-37,778,151", // is the deletion homozygous?
//        "--regions", "22:37,980,155-37,980,306", // false positive deletion
//        "--regions", "22:38,082,212-38,082,363", // false positive SNPs
//        "--regions", "22:38,092,378-38,092,418", // not left aligning the insertion at 38092390
//        "--regions", "22:39,517,648-39,517,952", // is the deletion & insetion homozygous?
//        "--regions", "22:40,946,793-40,946,944", // where is the deletion?
//        "--regions", "22:40,952,705-40,953,009", // another false positive deletion
//        "--regions", "22:41,434,235-41,434,310", // false positives
//        "--regions", "22:41,509,085-41,509,236", // false negative insertion (GATK calls, Platypus doesn't)?
//        "--regions", "22:41,988,173-41,988,324", // GATK is calling a deletion at 41988249...
        
        //"--regions", "22",
        
        //"--regions", "22:16,909,216-16,909,255", // NA12878_low_coverage insertion followed by 2 SNPs
        
        // Bad VCF representation
        //"--regions", "22:41,836,102-41,836,142", // NA12878_low_coverage deletion 2|0
        //"--regions", "22:43,035,829-43,035,869", // NA12878_low_coverage deletion 2|0
        
        //"--regions", "16:46,392,879-46,393,098", // NA12878_low_coverage huge memory spike
        //"--regions", "22:20656122-20656146", // NA12878_high_coverage causing memory spike
        //"--regions", "1:224,024,837-224,024,877", // NA12878_low_coverage edge case insertions
        //"--regions", "1:874950-874951", // NA12878_low_coverage edge case insertions
        //"--regions", "22:16,232,038-16,232,117", // NA12878_high_coverage VERY ineresting deletions!!
        //"--regions", "Y:13451184-13451197", // crazy Y haplotypes NA12878
        //"--regions", "22:37,276,961-37,277,045", // interesting dual deletions NA12878
        //"--regions", "22:23,474,658-23,475,205", // complex insertion/SNP NA12878
        //"--regions", "22:16,231,876-16,232,206", // Complex phasing region
        //"--regions", "22:16,139,705-16,140,266", // NA12878_high_coverage complex region
        //"-L", "6:31,235,688-31,235,883", // TEST phaser
        //"-L", "5:2,726,554-2,726,594", // snp just before deletion
        //"-L", "5:92,593,056", // SNP next to another SNP
        //"-L", "5:92,593,057", // SNP next to another SNP
        //"-L", "5:64,525,420-64,525,459", // two close SNPS
        
        //"--regions", "5:64,525,965-64,526,830", // NA12878_high_coverage crazy alignments
        
        //"--regions", "6:93,705,800-93,706,166", // NA12878_low_coverage no phase
        //"--regions", "6:58,605,652-58,606,061",   // NA12878_low_coverage partial phase
        //"--regions", "6:58,605,687-58,605,779",   // NA12878_low_coverage phase strong
        //"--regions", "3:108,803,741-108,803,854", // NA12878_low_coverage phase weak
        //"--regions", "6:22,877,929-22,878,012", // NA12878_low_coverage HMM error
        //"--regions", "6:89,236,310-89,237,082", // NA12878_high_coverage very nice phasing test
        //"--regions", "21:11,063,185-11,063,327", // weird haplotpes
        //"--regions", "3:47,251,793-47,251,839", // NA12878_high_coverage interesting indels
        //"--regions", "Y:13447283-13447483", // NA12878_low_coverage too many haplotypes
        
        // for cancer
        //"--regions", "6:52,873,970-52,882,816",
        //"--regions", "5:76,747,066-76,747,106", // not a reversion
        //"--regions", "5:76,781,703-76,781,743", // not a reversion
        //"--regions", "5:76,785,333-76,785,478", // HMM error
        
        //"--skip-regions", "1:1,000,000-2,000,000", "1:1,500,000-10,000,000",
        //"--skip-regions-file", human_skip_regions.c_str(),
        
        // candidate parameters
        "--min-supporting-reads", "2",
        "--min-snp-base-quality", "20",
        
        // read filters
        "--min-mapping-quality", "20",
        //"--allow-marked-duplicates",
        //"--allow-octopus-duplicates",
        "--trim-soft-clipped",
        //"--tail-trim-size", "3",
        "--trim-adapters",
        "--no-unmapped-mates",
        
        "--reference-cache-size", "20000",
        //"--no-downsampling",
        "--downsample-above", "300",
        "--downsample-target", "200",
        
        //"--no-candidates-from-alignments",
        //"--candidates-from-assembler",
        //"--kmer-size", "5",
        //"--candidates-from-source", sample_vcf.c_str(),
        //"--regenotype",
        
        "--min-variant-posterior", "10",
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
