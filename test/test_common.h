//
//  test_common.h
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_test_common_h
#define Octopus_test_common_h

namespace details {

    static std::string home_dir {getenv("HOME")};

    static std::string genomics_dir {"/Genomics/"};
    static std::string reference_dir {"References/"};
    static std::string bam_dir {"Illumina/"};

    // Reference names
    static std::string ecoli_reference_name {"R00000042.fasta"};
    static std::string human_reference_name {"human_g1k_v37.fasta"};
    static std::string lambda_reference_name {"lambda_ref.fasta"};

    // BAMs
    static std::string human_1000g_bam_name0 {"HG00100.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name1 {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name2 {"HG00102.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name3 {"HG00103.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"};
    static std::string ecoli_bam_name {"WTCHG_119208_201103.bam"};

    // CRAMs
    static std::string human_1000g_cram_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.cram"};
    
    // VCFs
    static std::string sample_vcf_name {"CEU.low_coverage.2010_07.xchr.sites.vcf.gz"};
    
    // donna data
    static std::string donna_dir {"donna/"};
    static std::string donna_bam_name1 {"356_005_sorted_with_labels.bam"};
    static std::string donna_bam_name2 {"356_006_sorted_with_labels.bam"};
    static std::string donna_bam_name3 {"357_005_sorted_with_labels.bam"};
    static std::string donna_bam_name4 {"357_006_sorted_with_labels.bam"};
}

// Full paths
static std::string human_reference_fasta {details::home_dir + details::genomics_dir +
                details::reference_dir + details::human_reference_name};

static std::string human_reference_fasta_index {human_reference_fasta + ".fai"};

static std::string ecoli_reference_fasta {details::home_dir + details::genomics_dir +
        details::reference_dir + details::ecoli_reference_name};

static std::string lambda_reference_fasta {details::home_dir + details::genomics_dir +
        details::reference_dir + details::lambda_reference_name};

static std::string human_1000g_bam0 {details::home_dir + details::genomics_dir + details::bam_dir +
    details::human_1000g_bam_name0};

static std::string human_1000g_bam1 {details::home_dir + details::genomics_dir + details::bam_dir +
        details::human_1000g_bam_name1};

static std::string human_1000g_bam2 {details::home_dir + details::genomics_dir + details::bam_dir +
        details::human_1000g_bam_name2};

static std::string human_1000g_bam3 {details::home_dir + details::genomics_dir + details::bam_dir +
    details::human_1000g_bam_name3};

static std::string human_1000g_cram {details::home_dir + details::genomics_dir + details::bam_dir +
        details::human_1000g_cram_name};

static std::string ecoli_bam {details::home_dir + details::genomics_dir + details::bam_dir +
    details::ecoli_bam_name};

static std::string sample_vcf {details::home_dir + details::genomics_dir + details::sample_vcf_name};

static std::string donna_bam1 {details::home_dir + details::genomics_dir + details::donna_dir +
    details::donna_bam_name1};
static std::string donna_bam2 {details::home_dir + details::genomics_dir + details::donna_dir +
    details::donna_bam_name2};
static std::string donna_bam3 {details::home_dir + details::genomics_dir + details::donna_dir +
    details::donna_bam_name3};
static std::string donna_bam4 {details::home_dir + details::genomics_dir + details::donna_dir +
    details::donna_bam_name4};

#endif
