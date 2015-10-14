//
//  test_common.h
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_test_common_h
#define Octopus_test_common_h

#include <vector>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace detail
{
    static std::string home_dir {getenv("HOME")};
    
    static std::string genomics_dir {"/Genomics/"};
    static std::string reference_dir {"References/"};
    static std::string bam_dir {"Illumina/"};
    static std::string octopus_test_dir {"/Development/Octopus/test/"};
    static std::string sample_vcf_dir {"sample_vcf/"};
    static std::string cancer_test_dir {"cancer/simulated/"};
    
    // Reference
    static std::string ecoli_reference_name {"R00000042.fasta"};
    static std::string human_reference_name {"human_g1k_v37.fasta"};
    static std::string lambda_reference_name {"lambda_ref.fasta"};
    
    // BAM
    
    // 1000G
    
    static std::string NA12878_low_coverage_name {"NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"};
    static std::string NA12878_high_coverage_name {"NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam"};
    static std::string NA12891_high_coverage_name {"NA12891.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam"};
    static std::string HG00101_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string HG00102_name {"HG00102.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    static std::string HG00103_name {"HG00103.mapped.ILLUMINA.bwa.GBR.low_coverage.20120522.bam"};
    
    // other
    
    static std::string ecoli_bam_name {"WTCHG_119208_201103.bam"};
    
    static std::string NA12878_simulated_cancer_basic {"NA12878_simulated_cancer_basic.bam"};
    
    // CRAM
    static std::string HG00101_cram_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.cram"};
    
    // VCF/BCF
    //static std::string sample_vcf_name {"test_triploid.vcf.gz"};
    static std::string sample_vcf_name {"CEU.low_coverage.2010_07.xchr.genotypes.vcf.gz"};
    static std::string sample_tabix_vcf_name {"CHBJPT.low_coverage.2010_07.xchr.sites.vcf.gz"};
    static std::string sample_bcf_name {"CHBJPT.low_coverage.2010_07.xchr.sites.vcf.gz"};
    //static std::string sample_vcf_name {"platypus.vcf"};
    
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
} // namespace detail

// Full paths

// regions

static const fs::path regions_txt_file {detail::home_dir + detail::octopus_test_dir + "test_regions.txt"};
static const fs::path regions_bed_file {detail::home_dir + detail::octopus_test_dir + "test_regions.bed"};
static const fs::path reads_file {detail::home_dir + detail::octopus_test_dir + "test_files.txt"};

// references

static const fs::path human_reference_fasta {detail::home_dir + detail::genomics_dir +
                    detail::reference_dir + detail::human_reference_name};

static const fs::path human_reference_fasta_index {human_reference_fasta.string() + ".fai"};

static const fs::path ecoli_reference_fasta {detail::home_dir + detail::genomics_dir +
        detail::reference_dir + detail::ecoli_reference_name};

static const fs::path lambda_reference_fasta {detail::home_dir + detail::genomics_dir +
        detail::reference_dir + detail::lambda_reference_name};

// reads

static const fs::path NA12878_low_coverage {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::NA12878_low_coverage_name};
    
static const fs::path NA12878_high_coverage {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::NA12878_high_coverage_name};

static const fs::path NA12891_high_coverage {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::NA12891_high_coverage_name};

static const fs::path HG00101 {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::HG00101_name};

static const fs::path HG00102 {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::HG00102_name};

static const fs::path HG00103 {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::HG00103_name};

static const fs::path HG00101_cram {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::HG00101_cram_name};

static const fs::path ecoli_bam {detail::home_dir + detail::genomics_dir + detail::bam_dir +
    detail::ecoli_bam_name};

static const std::vector<fs::path> bams_1000G {NA12878_low_coverage, };

// vcfs

static const fs::path sample_vcf {detail::home_dir + detail::genomics_dir  + detail::sample_vcf_dir +
    detail::sample_vcf_name};

static const fs::path test_out_vcf {detail::home_dir + detail::octopus_test_dir + "test.vcf"};
static const fs::path test_out_vcfgz {detail::home_dir + detail::octopus_test_dir + "test.vcf.gz"};
static const fs::path test_out_bcf {detail::home_dir + detail::octopus_test_dir + "test.bcf"};

// hiv

static const fs::path hiv_reference {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reference_dir + detail::hiv_reference_name};

static const fs::path hiv_bam1 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name1};
static const fs::path hiv_bam2 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name2};
static const fs::path hiv_bam3 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name3};
static const fs::path hiv_bam4 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name4};
static const fs::path hiv_bam5 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name5};
static const fs::path hiv_bam6 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name6};
static const fs::path hiv_bam7 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name7};
static const fs::path hiv_bam8 {detail::home_dir + detail::genomics_dir + detail::hiv_dir +
    detail::hiv_reads_dir + detail::hiv_bam_name8};

static const std::vector<fs::path> hiv_bams {hiv_bam1, hiv_bam2, hiv_bam3, hiv_bam4, hiv_bam5, hiv_bam6, hiv_bam7, hiv_bam8};

// simulated cancer

static const fs::path NA12878_simulated_cancer_basic {detail::home_dir + detail::genomics_dir +
    detail::cancer_test_dir + "NA12878.13.36802245-36805221.final.bam"};

#endif
