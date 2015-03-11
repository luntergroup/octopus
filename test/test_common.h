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
    static std::string human_1000g_bam_name1 {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name2 {"HG00472.mapped.ILLUMINA.bwa.CHS.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name3 {"HG00705.mapped.ILLUMINA.bwa.CHS.low_coverage.20120522.bam"};
    static std::string chrom_20_bam_name {"human_chr_20.bam"};

    // CRAMS
    static std::string human_1000g_cram_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.cram"};
        
}

// Full paths
static std::string human_reference_fasta {details::home_dir + details::genomics_dir +
                details::reference_dir + details::human_reference_name};

static std::string human_reference_fasta_index {human_reference_fasta + ".fai"};

static std::string ecoli_reference_fasta {details::home_dir + details::genomics_dir +
        details::reference_dir + details::ecoli_reference_name};

static std::string lambda_reference_fasta {details::home_dir + details::genomics_dir +
        details::reference_dir + details::lambda_reference_name};

static std::string human_1000g_bam1 {details::home_dir + details::genomics_dir + details::bam_dir +
        details::human_1000g_bam_name1};

static std::string human_1000g_bam2 {details::home_dir + details::genomics_dir + details::bam_dir +
        details::human_1000g_bam_name2};

static std::string human_1000g_bam3 {details::home_dir + details::genomics_dir + details::bam_dir +
    details::human_1000g_bam_name3};

static std::string human_1000g_cram {details::home_dir + details::genomics_dir + details::bam_dir +
        details::human_1000g_cram_name};

static std::string human_1000g_bam1_chrom_20 {details::home_dir + details::genomics_dir + details::bam_dir +
        details::chrom_20_bam_name};

#endif
