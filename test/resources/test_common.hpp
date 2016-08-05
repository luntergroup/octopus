// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef Octopus_test_common_h
#define Octopus_test_common_h

#include <vector>
#include <boost/filesystem.hpp>

namespace fs = boost::filesystem;

namespace octopus { namespace test { namespace detail {
    // common directories
    
    static std::string home_dir {getenv("HOME")};
    static std::string genomics_dir {"/Genomics/"};
    static std::string reference_dir {"References/"};
    static std::string bam_dir {"Illumina/"};
    static std::string octopus_test_dir {"/Development/octopus/test/"};
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
    
    // CRAM
    static std::string HG00101_cram_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.cram"};
    static std::string NA12878_low_coverage_cram_name {"NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.cram"};
    
    // VCF/BCF
    
    static std::string sample_vcf_name {"CEU.low_coverage.2010_07.xchr.genotypes.vcf"};
    static std::string sample_vcfgz_name {"CEU.low_coverage.2010_07.xchr.genotypes.vcf.gz"};
    static std::string sample_bcf_name {"CEU.low_coverage.2010_07.xchr.genotypes.bcf"};
} // namespace detail

// Full paths

// regions

static const fs::path regions_txt_file {detail::home_dir + detail::octopus_test_dir + "test_regions.txt"};
static const fs::path regions_bed_file {detail::home_dir + detail::octopus_test_dir + "test_regions.bed"};
static const fs::path reads_file {detail::home_dir + detail::octopus_test_dir + "test_files.txt"};
static const fs::path human_skip_regions {detail::home_dir + detail::octopus_test_dir + "human_skip_regions2.txt"};
static const fs::path hla_regions {detail::home_dir + detail::octopus_test_dir + "hla_regions.txt"};

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
static const fs::path NA12878_low_coverage_cram {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::NA12878_low_coverage_cram_name};
static const fs::path ecoli_bam {detail::home_dir + detail::genomics_dir + detail::bam_dir +
        detail::ecoli_bam_name};

// vcfs

static const fs::path sample_vcf {detail::home_dir + detail::genomics_dir  + detail::sample_vcf_dir +
    detail::sample_vcf_name};
static const fs::path sample_vcfgz {detail::home_dir + detail::genomics_dir  + detail::sample_vcf_dir +
    detail::sample_vcfgz_name};
static const fs::path sample_bcf {detail::home_dir + detail::genomics_dir  + detail::sample_vcf_dir +
    detail::sample_bcf_name};
static const fs::path test_out_vcf {detail::home_dir + detail::octopus_test_dir + "test.vcf"};
static const fs::path test_out_vcfgz {detail::home_dir + detail::octopus_test_dir + "test.vcf.gz"};
static const fs::path test_out_bcf {detail::home_dir + detail::octopus_test_dir + "test.bcf"};

// misc

static const fs::path non_existent {detail::home_dir + detail::octopus_test_dir + "non_existent"};
static const fs::path contig_ploidies_txt_file {detail::home_dir + detail::octopus_test_dir + "contig_ploidies.txt"};

// Helper methods

inline bool test_file_exists(const fs::path& file)
{
    return fs::exists(file);
}

inline bool all_region_test_files_exist()
{
    return test_file_exists(regions_txt_file) && test_file_exists(regions_bed_file)
            && test_file_exists(reads_file) && test_file_exists(human_skip_regions);
}

inline bool all_reference_test_files_exist()
{
    return test_file_exists(human_reference_fasta) && test_file_exists(ecoli_reference_fasta)
            && test_file_exists(lambda_reference_fasta);
}

inline bool all_bam_test_files_exist()
{
    return test_file_exists(NA12878_low_coverage) && test_file_exists(NA12878_high_coverage)
            && test_file_exists(NA12891_high_coverage);
}

inline bool all_cram_test_files_exist()
{
    return test_file_exists(NA12878_low_coverage_cram);
}

inline bool all_read_test_files_exist()
{
    return all_bam_test_files_exist() && all_cram_test_files_exist();
}

inline bool all_variant_test_files_exist()
{
    return test_file_exists(sample_vcf);
}

inline bool all_misc_files_exist()
{
    return test_file_exists(non_existent) && test_file_exists(contig_ploidies_txt_file);
}

inline bool all_test_files_exist()
{
    return     all_region_test_files_exist()
            && all_reference_test_files_exist()
            && all_read_test_files_exist()
            && all_variant_test_files_exist()
            && all_misc_files_exist();
}

inline bool is_writable_test_file(const fs::path& file)
{
    return file == test_out_vcf || file == test_out_vcfgz || file == test_out_bcf;
}

inline void remove_test_file(const fs::path& file)
{
    if (is_writable_test_file(file)) {
        fs::remove(file);
    }
}

inline void cleanup_test_files()
{
    remove_test_file(test_out_vcf);
    remove_test_file(test_out_vcfgz);
    remove_test_file(test_out_bcf);
}

} // namespace test
} // namespace octopus

#endif
