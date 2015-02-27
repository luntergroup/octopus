//
//  test_common.h
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_test_common_h
#define Octopus_test_common_h

#include "catch.hpp"
#define CATCH_CONFIG_MAIN

static std::string home_dir {getenv("HOME")};

static std::string genomics_dir {"/Genomics/"};
static std::string reference_dir {"References/"};
static std::string bam_dir {"Illumina/"};

// Reference names
static std::string ecoli_reference {"R00000042.fasta"};
static std::string human_reference {"human_g1k_v37.fasta"};
static std::string lambda_reference {"lambda_ref.fasta"};

// BAMs
static std::string human_1000g_bam_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};

// CRAMS
static std::string human_1000g_cram_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.cram"};

// Full paths
static std::string human_reference_fasta {home_dir + genomics_dir + reference_dir + human_reference};
static std::string human_reference_fasta_index {human_reference_fasta + ".fai"};
static std::string ecoli_reference_fasta {home_dir + genomics_dir + reference_dir + ecoli_reference};
static std::string lambda_reference_fasta {home_dir + genomics_dir + reference_dir + lambda_reference};

static std::string human_1000g_bam {home_dir + genomics_dir + bam_dir + human_1000g_bam_name};
static std::string human_1000g_cram {home_dir + genomics_dir + bam_dir + human_1000g_cram_name};

#endif
