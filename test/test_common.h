//
//  test_common.h
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_test_common_h
#define Octopus_test_common_h

namespace detail
{
    static std::string home_dir {getenv("HOME")};
    
    static std::string genomics_dir {"/Genomics/"};
    static std::string reference_dir {"References/"};
    static std::string bam_dir {"Illumina/"};
    static std::string octopus_test_dir {"/Development/Octopus/test/"};
    static std::string sample_vcf_dir {"sample_vcf/"};
    
    // Reference
    static std::string ecoli_reference_name {"R00000042.fasta"};
    static std::string human_reference_name {"human_g1k_v37.fasta"};
    static std::string lambda_reference_name {"lambda_ref.fasta"};
    
    // BAM
    static std::string human_1000g_bam_name0 {"HG00100.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name1 {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name2 {"HG00102.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string human_1000g_bam_name3 {"HG00103.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"};
    static std::string ecoli_bam_name {"WTCHG_119208_201103.bam"};
    
    // CRAM
    static std::string human_1000g_cram_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.cram"};
    
    // VCF/BCF
    //static std::string sample_vcf_name {"test_triploid.vcf.gz"};
    static std::string sample_vcf_name {"CEU.low_coverage.2010_07.xchr.genotypes.vcf.gz"};
    static std::string sample_tabix_vcf_name {"CHBJPT.low_coverage.2010_07.xchr.sites.vcf.gz"};
    static std::string sample_bcf_name {"CHBJPT.low_coverage.2010_07.xchr.sites.vcf.gz"};
    //static std::string sample_vcf_name {"platypus.vcf"};
    
    // donna data
    static std::string donna_dir {"donna/"};
    static std::string donna_bam_name1 {"356_005_sorted_with_labels.bam"};
    static std::string donna_bam_name2 {"356_006_sorted_with_labels.bam"};
    static std::string donna_bam_name3 {"357_005_sorted_with_labels.bam"};
    static std::string donna_bam_name4 {"357_006_sorted_with_labels.bam"};
    
    // hiv data
    static std::string hiv_dir {"hiv/"};
    static std::string hiv_reference_dir {"references/"};
    static std::string hiv_reads_dir {"reads/alignments/"};
    static std::string hiv_reference_name {"B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta"};
    static std::string hiv_bam_name1 {"10065_1_28_1_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
    static std::string hiv_bam_name2 {"10065_1_28_2_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
    static std::string hiv_bam_name3 {"12370_1_1_1_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
    static std::string hiv_bam_name4 {"12370_1_1_2_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
    static std::string hiv_bam_name5 {"12426_1_49_1_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
    static std::string hiv_bam_name6 {"12426_1_49_2_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
    static std::string hiv_bam_name7 {"13591_1_16_1_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
    static std::string hiv_bam_name8 {"13591_1_16_2_B.FR.83.HXB2_LAI_IIIB_BRU.K03455_bwa.bam"};
} // end namespace detail

// Full paths

// regions

static std::string regions_txt_file {detail::home_dir + detail::octopus_test_dir + "test_regions.txt"};
static std::string regions_bed_file {detail::home_dir + detail::octopus_test_dir + "test_regions.bed"};
static std::string reads_file {detail::home_dir + detail::octopus_test_dir + "test_files.txt"};

// references

static std::string human_reference_fasta {detail::home_dir + detail::genomics_dir +
                detail::reference_dir + detail::human_reference_name};

static std::string human_reference_fasta_index {human_reference_fasta + ".fai"};

static std::string ecoli_reference_fasta {detail::home_dir + detail::genomics_dir +
        detail::reference_dir + detail::ecoli_reference_name};

static std::string lambda_reference_fasta {detail::home_dir + detail::genomics_dir +
        detail::reference_dir + detail::lambda_reference_name};

// reads

static std::string human_1000g_bam0 {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::human_1000g_bam_name0};

static std::string human_1000g_bam1 {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::human_1000g_bam_name1};

static std::string human_1000g_bam2 {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::human_1000g_bam_name2};

static std::string human_1000g_bam3 {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::human_1000g_bam_name3};

static std::string human_1000g_cram {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::human_1000g_cram_name};

static std::string ecoli_bam {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::ecoli_bam_name};

// vcfs

static std::string sample_vcf {detail::home_dir + detail::genomics_dir  + detail::sample_vcf_dir +
    detail::sample_vcf_name};

// donna

static std::string donna_bam1 {detail::home_dir + detail::genomics_dir + detail::donna_dir +
    detail::donna_bam_name1};
static std::string donna_bam2 {detail::home_dir + detail::genomics_dir + detail::donna_dir +
    detail::donna_bam_name2};
static std::string donna_bam3 {detail::home_dir + detail::genomics_dir + detail::donna_dir +
    detail::donna_bam_name3};
static std::string donna_bam4 {detail::home_dir + detail::genomics_dir + detail::donna_dir +
    detail::donna_bam_name4};

// hiv

static std::string hiv_reference {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reference_dir + detail::hiv_reference_name};

static std::string hiv_bam1 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name1};
static std::string hiv_bam2 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name2};
static std::string hiv_bam3 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name3};
static std::string hiv_bam4 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name4};
static std::string hiv_bam5 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name5};
static std::string hiv_bam6 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name6};
static std::string hiv_bam7 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name7};
static std::string hiv_bam8 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name8};

#endif
