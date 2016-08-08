// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef Octopus_mock_options_h
#define Octopus_mock_options_h

#include <boost/program_options.hpp>
#include <boost/optional.hpp>

#include "option_parser.hpp"
#include "test_common.hpp"

namespace octopus { namespace test {

namespace po = boost::program_options;

inline auto get_basic_mock_options()
{
    const char *argv[] = {
        "octopus",
        
        //"--help",
        //"--version",
        
        "--debug",
        //"--trace",
        
        //"--sites-only",
        
        //"--samples", "NOT-A-SAMPLE",
        
        "--working-directory", "~/Genomics/giab_benchmark",
        //"--working-directory", "~/Genomics/MCG",
        //"--working-directory", "~/Genomics/cancer/TCGA/benchmark",
        
        //"--target-read-buffer-size", "1.0",
        //"--reference-cache-size", "100",
        //"--threads",
        
        //"--contig-output-order", "AsInReferenceReversed",
        
        "--reference", human_reference_fasta.c_str(),
        //"--reference", ecoli_reference_fasta.c_str(),
        
        "--reads", NA12878_low_coverage.c_str(),
        //"--reads", NA12878_high_coverage.c_str(),
        //"--reads", "~/Genomics/Illumina/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.chr22.bam",
        
        //"--reads", NA12891_high_coverage.c_str(),
        //"--reads", HG00101.c_str(),
        
        //"--reads", NA12878_low_coverage.c_str(), HG00101.c_str(), HG00102.c_str(), HG00103.c_str(),
        //"--reads", NA12878_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12891_high_coverage.c_str(), NA12878_simulated_cancer_basic.c_str(), //cancer test
        //"--reads", NA12878_high_coverage.c_str(), NA12891_high_coverage.c_str(),
        
        //"--reads", ecoli_bam.c_str(),
        
        //"--reads", "~/Genomics/Illumina/NA12878HC_HLA_C.bam",
        
        // TCGA
        //"--reads", "~/Genomics/cancer/TCGA/benchmark/G15511.HCC1143_BL.1.chr22.bam",
        //"--reads", "~/Genomics/cancer/TCGA/benchmark/G15511.HCC1143.1.chr22.bam",
        
        // MCG
        //"--reads", "~/Genomics/MCG/10120_chr2_47641558_GTA_G.RG.bam",
        //"--reads", "~/Genomics/MCG/D59597_Cov3.RG.bam",
        
        //"--caller", "cancer", // default "population"
        //"--normal-sample", "HCC1143 BL",
        
        //"--organism-ploidy", "3",
        //"--contig-ploidies", "MT=1", "Y=1",
        //"--contig-ploidies-file", contig_ploidies_txt_file.c_str(),
        
        //"--make-blocked-refcalls",
        //"--make-positional-refcalls",
        
        /* input regions */
        
        //"--use-one-based-indexing",
        
        "--regions", "22:30,000,000-30,010,000",
        
        //"--regions", "22",
        
        //"--regions", "2:120,098,309-120,098,368", // NA12878 LC - GQ is too high
        
        //"--regions", "22:22,583,864-22,584,168", // FN
        //"--regions", "22:40,015,775-40,016,384", // FP
        //"--regions", "22:48,744,696-48,744,879", // difficult region (FPs)
        //"--regions", "22:47,416,419-47,417,744", // bad model filtering
        
        //"--regions-file", "data/NA12878_GIAB_highconf_regions_big_expanded.bed",
        //"--regions-file", "data/NA12878_GIAB_chr22_highconf_regions_big_expanded.bed",
        
        //"--regions", "22:19,659,923-19,661,143", // very cool assembler deletion
        //"--regions", "22:49,217,985-49,218,594", // very cool assembler deletion
        //"--regions", "22:24,712,696-24,713,000", // assembler proposing insane deletion
        
        //"--regions", "22:48744706-48744869", // tricky GIAB hcr
        
        //"--regions", "22:44,695,012-44,695,163", // alignment routine going over bounds causing FP
        //"--regions", "22:35,392,836-35,393,140", // alignment routine going over bounds causing wrong genotype call
        //"--regions", "22:33,862,663-33,862,749", // SNV error model causing FN
        //"--regions", "22:47,711,210-47,711,361", // SNV error model causing low qual TP
        
        //"--regions", "22:42,943,057-42,943,471", // why is GATK calling these?
        
        //"--regions", "6:31,236,339-31,240,057", // HLA-C
        
        // MCG
        //"--regions", "2:47,640,494-47,644,304", // whole region
        //"--regions", "2:47,643,156-47,643,222", // true snp
        //"--regions", "2:47,640,572-47,642,302", // first block (no variants)
        //"--regions", "2:47,641,997-47,644,890", // second block (one SNP)
        
        // for cancer
        
        // TCGA HCC1143
        
        //"--regions", "22:29,606,605-29,606,972", // Bad model filter!
        
        //"--regions", "22:18,444,570-18,444,787", // real somatic SNV?
        //"--regions", "22:25,340,394-25,340,434", // real somatic SNV?
        //"--regions", "22:28,681,165-28,681,395", // real somatic SNV?
        //"--regions", "22:27,707,987-27,708,223", // FN low frequency somatic deletion?
        //"--regions", "22:37,166,974-37,167,210", // FN real somatic deletion?
        //"--regions", "22:28,553,936-28,554,100", // somatic insertion? Low posterior
        //"--regions", "22:28,645,596-28,645,926", // Real somatic SNV? Low posterior
        //"--regions", "22:28,970,260-28,970,424", // real somatic deletion?
        //"--regions", "22:29,474,547-29,474,587", // real somatic SNV?
        //"--regions", "22:29,761,141-29,761,305", // real somatic SNV?
        
        //"--regions", "22:46,729,130-46,730,350", // PASSing FP SNV
        //"--regions", "22:40,946,818-40,946,893", // PASSing FP SNV
        //"--regions", "22:41,607,540-41,607,580", // tricky PASSing FP SNVs
        //"--regions", "22:42,905,833-42,905,873", // PASSing FP SNV
        //"--regions", "22:43,410,348-43,410,423", // PASSing FP SNV
        //"--regions", "22:43,622,482-43,622,832", // PASSing FP SNV. Fixed: haplotype filtering
        //"--regions", "22:47,094,596-47,094,671", // PASSing FP SNV
        //"--regions", "22:47,996,352-47,996,392", // PASSing FP SNV
        
        //"--regions", "22:40,000,000-41,000,000",
        
//        "--regions", "22:41,120,796-41,120,871", // FP SNP. Strand bias. Run through.
//        "--regions", "22:43,960,301-43,960,376", // FP SNP. Strand bias. Run through.
//        "--regions", "22:42,508,653-42,508,957", // FP SNP. Strand bias. Low MQ.
//        "--regions", "22:42,726,975-42,727,015", // FP SNP. Strand bias.
//        "--regions", "22:43,468,078-43,468,229", // FP SNPs. Strand bias. Run through.
//        "--regions", "22:44,202,032-44,202,107", // FP SNPs. Strand bias. Run through.
//        "--regions", "22:44,230,824-44,230,975", // FP SNP. MQ. Strand bias.
//        "--regions", "22:44,243,370-44,243,573", // FP SNP. Strand bias. Run through.
//        "--regions", "22:45,332,712-45,332,787", // FP SNP. Strand bias. Run through.
//        "--regions", "22:45,895,589-45,895,664", // FP SNP. Strand bias. Run through.
//        "--regions", "22:46,364,941-46,365,092", // FP SNP. Strand bias.
//        "--regions", "22:46,491,358-46,491,509", // FP SNP. Low MQ.
//        "--regions", "22:49,872,443-49,872,607", // FP SNP. Strand bias.
//        "--regions", "22:46,729,658-46,729,822", // FP SNP. Strand bias
//        "--regions", "22:46,491,424-46,491,728", // FP SNP
//        "--regions", "22:46,614,813-46,614,853", // FP SNP. Strand bias.
//        "--regions", "22:46,729,664-46,729,815", // FP SNP. Strand bias.
        
        //"--regions", "22:43,445,777-43,445,817", // FP del
        //"--regions", "22:42,937,177-42,937,252", // FP ins
        //"--regions", "22:40,603,966-40,604,296", // FP del
        //"--regions", "22:30,213,959-30,214,123", // FP somatic SNV. Very difficult.
        //"--regions", "22:33,595,824-33,595,864", // FP somatic del
        //"--regions", "22:25,731,146-25,731,476", // FP SNP. What's going on here?
        //"--regions", "22:28,063,906-28,064,070", // What's going on here?
        
        //"--regions", "22:27,297,910-27,297,950", // Not a FP, somaticSniper is calling, but what is going on?
        
        //"--regions", "22:41,528,865-41,529,474", // FP ins. Solved: Conservative phasing
        //"--regions", "22:45,748,970-45,749,325", // FP del. Solved: up max-haplotypes (true haplotype getting filtered)
        //"--regions", "22:25,656,657-25,656,697", // FP SNP. Solved: Conservative phasing
        //"--regions", "22:43,358,707-43,358,858", // FP ins. Solved: Conservative phasing
        //"--regions", "22:42,508,785-42,508,825", // FP SNP. Solved: Disable MQ filter
        //"--regions", "22:28,812,452-28,812,616", // FP del. Solved: indel error model
        //"--regions", "22:29,455,729-29,455,769", // FP del. Solved: indel error model
        //"--regions", "22:23,842,194-23,842,524", // FP SNV. Solved: Disable overlap masking. Why?
        //"--regions", "22:23,869,889-23,870,053", // FP del. Solved: Disable marked dups filter
        //"--regions", "22:24,160,300-24,160,359", // FP del. Solved: Disable soft clipping
        //"--regions", "22:25,055,737-25,055,901", // FP del. Solved: MQ Filter
        //"--regions", "22:25,656,651-25,656,728", // FP SNP. Solved: Disable MQ filter or allow lagging
        //"--regions", "22:25,330,405-25,330,569", // FP DEL. Solved: model filter
        //"--regions", "22:25,070,597-25,070,854", // FP SNP. Solved: model filter
        
        // For ecoli
        //"--regions", "R00000042:3008660-3020000",
        
        //"--regions-file", hla_regions.c_str(),
        
        //"--skip-regions", "1:1,000,000-2,000,000", "1:1,500,000-10,000,000",
        //"--skip-regions-file", human_skip_regions.c_str(),
        
        //"--regenotype", "~/Genomics/octopus_test/AllVariants.vcf",
        //"--regenotype", "~/Genomics/cancer/TCGA/benchmark/somatic_sniper_chr22_Q40.vcf",
        
        // read transforms
        //"--disable-read-transforms",
        //"--disable-soft-clip-masking",
        //"--mask-soft-clipped-boundries", "0",
        //"--disable-adapter-masking",
        //"--disable-overlap-masking",
        
        // read filters
        //"--disable-read-filtering",
        //"--consider-unmapped-reads",
        //"--consider-reads-with-unmapped-segments",
        //"--min-mapping-quality", "1",
        //"--good-base-quality", "0",
        //"--min-good-bases", "0",
        //"--allow-marked-duplicates",
        //"--allow-octopus-duplicates",
        //"--allow-qc-fails",
        //"--consider-reads-with-unmapped-segments",
        //"--allow-adapter-contaminated-reads",
        
        // downsampling
        //"--disable-downsampling",
        //"--downsample-above", "300",
        //"--downsample-target", "200",
        
        // candidate generation
        //"--min-supporting-reads=1",
        //"--min-base-quality", "15",
        //"--disable-raw-cigar-candidate-generator",
        //"--kmer-size", "5",
        //"--candidates-from-source", sample_vcf.c_str(),
        //"--disable-assembly-candidate-generator",
        //"--kmer-size", "75",
        //"--assembler-mask-base-quality=2",
        //"--max-variant-size", "25",
        
        //"--phasing-level", "Minimal", // Minimal, Conservative, Aggressive
        //"--disable-inactive-flank-scoring",
        //"--max-haplotypes", "50",
        //"--min-haplotype-posterior", "1e-30",
        
        // filtering
        //"--disable-model-filtering",
        "--disable-call-filtering",
        
        //"--min-variant-posterior", "0.5",
        //"--min-refcall-posterior", "0",
        
        //"--snp-heterozygosity", "0.000333",
        //"--indel-heterozygosity", "0.01",
        
        // cancer
        //"--somatic-mutation-rate", "0.001",
        //"--min-somatic-frequency", "0.01",
        //"--credible-mass", "0.99",
        //"--min-somatic-posterior", "2",
        
        //"--somatics-only",
        
        //"--max-open-read-files", "1",
        
        //"--output", "octopus_hla.vcf",
        //"--output", "octopus_mcg.vcf",
        //"--output", "octopus_cancer.vcf",
        "--output", "octopus_calls_debug.vcf",
        //"--output", "octopus_calls_assemble.vcf",
        //"--output", "octopus_calls_fast.vcf",
        //"--output", "octopus/octopus_NA12878_LC.vcf.gz",
        
        //"--legacy",
        
        nullptr
    };
    
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    return octopus::options::parse_options(argc, argv);
}

} // namespace test
} // namespace octopus

#endif
