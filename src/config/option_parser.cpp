// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "option_parser.hpp"

#include <vector>
#include <array>
#include <regex>
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <fstream>

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/lexical_cast.hpp>

#include "utils/path_utils.hpp"
#include "utils/memory_footprint.hpp"
#include "utils/string_utils.hpp"
#include "basics/phred.hpp"
#include "exceptions/user_error.hpp"
#include "config.hpp"

#include "core/csr/measures/measure_factory.hpp"

namespace po = boost::program_options;
namespace fs = boost::filesystem;

namespace octopus { namespace options {

fs::path resolve_path(const fs::path& path, const OptionMap& options);
void parse_config_file(const fs::path& config_file, OptionMap& vm, const po::options_description& options);

// boost::option cannot handle option dependencies so we must do our own checks
void conflicting_options(const OptionMap& vm, const std::string& opt1, const std::string& opt2);
void option_dependency(const OptionMap& vm, const std::string& for_what, const std::string& required_option);
void check_positive(const std::string& option, const OptionMap& vm);
void check_reads_present(const OptionMap& vm);
void check_region_files_consistent(const OptionMap& vm);
void check_trio_consistent(const OptionMap& vm);
void validate_caller(const OptionMap& vm);
void validate(const OptionMap& vm);

po::parsed_options run(po::command_line_parser& parser);

OptionMap parse_options(const int argc, const char** argv)
{
    po::options_description general("General");
    general.add_options()
    ("help,h",
     "Report detailed option information")
    
    ("version",
     "Report detailed version information")
    
    ("config",
     po::value<fs::path>(),
     "Config file to populate command line options")
    
    ("debug",
     po::value<fs::path>()->implicit_value("octopus_debug.log"),
     "Create log file for debugging")
    
    ("trace",
     po::value<fs::path>()->implicit_value("octopus_trace.log"),
     "Create very verbose log file for debugging")
    
    ("working-directory,w",
     po::value<fs::path>(),
     "Sets the working directory")
    
    ("resolve-symlinks",
     po::bool_switch()->default_value(false),
     "Replace all symlinks to their resolved targets")
    
    ("threads,tentacles",
     po::value<int>()->implicit_value(0),
     "Maximum number of threads to be used. If no argument is provided unlimited threads are assumed")
    
    ("max-reference-cache-memory,X",
     po::value<MemoryFootprint>()->default_value(*parse_footprint("500MB"), "500MB"),
     "Maximum memory for cached reference sequence")
    
    ("target-read-buffer-memory,B",
     po::value<MemoryFootprint>()->default_value(*parse_footprint("6GB"), "6GB"),
     "None-binding request to limit the memory of buffered read data")
    
    ("target-working-memory",
     po::value<MemoryFootprint>(),
     "Target working memory per thread for computation, not including read or reference data")
     
    ("max-open-read-files",
     po::value<int>()->default_value(250),
     "Limits the number of read files that are open simultaneously")

    ("temp-directory-prefix",
     po::value<fs::path>()->default_value("octopus-temp"),
     "File name prefix of temporary directory for calling")
    
    ("keep-temporary-files",
     po::bool_switch()->default_value(false),
     "Do not remove temporary files, even after an error")
    
    ("reference,R",
     po::value<fs::path>()->required(),
     "Indexed FASTA format reference genome file to be analysed")
    
    ("reads,I",
     po::value<std::vector<fs::path>>()->multitoken(),
     "Indexed BAM/CRAM files to be analysed")
    
    ("reads-file,i",
     po::value<std::vector<fs::path>>()->multitoken(),
     "Files containing lists of BAM/CRAM files, one per line, to be analysed")
    
    ("regions,T",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-separated list of regions (chrom:begin-end) to be analysed")
    
    ("regions-file,t",
     po::value<fs::path>(),
     "File containing a list of regions (chrom:begin-end), one per line, to be analysed")
    
    ("skip-regions,K",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-separated list of regions (chrom:begin-end) to skip")
    
    ("skip-regions-file,k",
     po::value<fs::path>(),
     "File of regions (chrom:begin-end), one per line, to skip")

    ("one-based-indexing",
     po::bool_switch()->default_value(false),
     "Assume one-based indexing rather than zero-based for input region options")
    
    ("samples,S",
     po::value<std::vector<std::string>>()->multitoken(),
     "Sub-set of sample names present in input read files to analyse")
    
    ("samples-file,s",
     po::value<fs::path>(),
     "File of sample names to analyse, one per line")
    
    ("ignore-unmapped-contigs",
     po::bool_switch()->default_value(false),
     "Ignore any contigs that are not mapped in the read files")
    
    ("pedigree",
     po::value<fs::path>(),
     "PED file containing sample pedigree")
    
    ("output,o",
     po::value<fs::path>(),
     "File to where output is written (calls are written to stdout if unspecified)")
    
    ("contig-output-order",
     po::value<ContigOutputOrder>()->default_value(ContigOutputOrder::referenceIndex),
     "The order that contigs should be written to the output [LEXICOGRAPHICAL_ASCENDING, LEXICOGRAPHICAL_DESCENDING, CONTIG_SIZE_ASCENDING, CONTIG_SIZE_DESCENDING, REFERENCE_INDEX, REFERENCE_INDEX_REVERSED]")
    
    ("sites-only",
     po::bool_switch()->default_value(false),
     "Only reports call sites (i.e. drop sample genotype information)")
    
    // ("regenotype",
    //  po::value<fs::path>(),
    //  "VCF file specifying calls to regenotype, only sites in this files will appear in the"
    //  " final output")
    
    ("bamout",
     po::value<fs::path>(),
     "Output realigned BAM files")
    
    ("bamout-type",
     po::value<RealignedBAMType>()->default_value(RealignedBAMType::mini),
     "Type of realigned evidence BAM to output [MINI, FULL]")
     
    ("data-profile",
     po::value<fs::path>(),
     "Output a profile of variation and errors found in the data")
    
    ("fast",
     po::bool_switch()->default_value(false),
     "Turns off some features to improve runtime, at the cost of worse calling accuracy and phasing")
    
    ("very-fast",
     po::bool_switch()->default_value(false),
     "Like --fast but even faster")
    ;
    
    po::options_description read_preprocessing("Read preprocessing");
    read_preprocessing.add_options()
    ("disable-read-preprocessing",
     po::bool_switch()->default_value(false),
     "Disable all read preprocessing")
    
    ("max-base-quality",
     po::value<int>(),
     "Cap all base qualities to this value")
     
    ("mask-low-quality-tails",
     po::value<int>()->implicit_value(3),
     "Masks read tail bases with base quality less than this")
     
    ("mask-tails",
     po::value<int>()->implicit_value(1),
     "Unconditionally mask this many read tail bases")
     
    ("mask-soft-clipped-bases",
     po::bool_switch()->default_value(false),
     "Enable masking of soft clipped bases")
    
    ("soft-clip-mask-threshold",
     po::value<int>()->implicit_value(3),
     "Only soft clipped bases with quality less than this will be recalibrated, rather than all bases")
    
    ("mask-soft-clipped-boundary-bases",
     po::value<int>()->default_value(2),
     "Masks this number of adjacent non soft clipped bases when soft clipped bases are present")
    
    ("disable-adapter-masking",
     po::bool_switch()->default_value(false),
     "Disable adapter detection and masking")
    
    ("disable-overlap-masking",
     po::bool_switch()->default_value(false),
     "Disable read segment overlap masking")
     
    ("mask-inverted-soft-clipping",
     po::bool_switch()->default_value(false),
     "Mask soft clipped sequence that is an inverted copy of a proximate sequence")
    
    ("mask-3prime-shifted-soft-clipped-heads",
     po::bool_switch()->default_value(false),
     "Mask soft clipped read head sequence that is a copy of a proximate 3' sequence")

    ("split-long-reads",
     po::bool_switch()->default_value(false),
     "Split reads longer than 'max-read-length' into linked fragments")
    
    ("consider-unmapped-reads",
     po::bool_switch()->default_value(false),
     "Allows reads marked as unmapped to be used for calling")
    
    ("min-mapping-quality",
     po::value<int>()->default_value(5),
     "Minimum read mapping quality required to consider a read for calling")
    
    ("good-base-quality",
     po::value<int>()->default_value(20),
     "Base quality threshold used by min-good-bases and min-good-base-fraction filters")
    
    ("min-good-base-fraction",
     po::value<double>()->implicit_value(0.5),
     "Base quality threshold used by min-good-bases filter")
    
    ("min-good-bases",
     po::value<int>()->default_value(20),
     "Minimum number of bases with quality min-base-quality before read is considered")
    
    ("allow-qc-fails",
     po::bool_switch()->default_value(false),
     "Filters reads marked as QC failed")
    
    ("min-read-length",
     po::value<int>(),
     "Filters reads shorter than this")
    
    ("max-read-length",
     po::value<int>()->default_value(10'000),
     "Filter reads longer than this")
    
    ("allow-marked-duplicates",
     po::bool_switch()->default_value(false),
     "Allows reads marked as duplicate in alignment record")
    
    ("allow-octopus-duplicates",
     po::bool_switch()->default_value(false),
     "Allows reads considered duplicates by octopus")
    
    ("duplicate-read-detection-policy",
     po::value<ReadDeduplicationDetectionPolicy>()->default_value(ReadDeduplicationDetectionPolicy::relaxed),
     "Policy to use for duplicate read detection [RELAXED, AGGRESSIVE]")
    
    ("allow-secondary-alignments",
     po::bool_switch()->default_value(false),
     "Allows reads marked as secondary alignments")
    
    ("allow-supplementary-alignments",
     po::bool_switch()->default_value(false),
     "Allows reads marked as supplementary alignments")
    
    ("no-reads-with-unmapped-segments",
     po::bool_switch()->default_value(false),
     "Filter reads with unmapped template segments to be used for calling")
    
    ("no-reads-with-distant-segments",
     po::bool_switch()->default_value(false),
     "Filter reads with template segments that are on different contigs")
    
    ("no-adapter-contaminated-reads",
     po::bool_switch()->default_value(false),
     "Filter reads with possible adapter contamination")
    
    ("max-decoy-supplementary-alignment-mapping-quality",
     po::value<int>()->default_value(5),
     "Filter reads with supplementary alignments mapped to decoy contigs with mapping quality greater than this")
    
    ("max-unplaced-supplementary-alignment-mapping-quality",
     po::value<int>()->default_value(5),
     "Filter reads with supplementary alignments mapped to unplaced contigs with mapping quality greater than this")
    
    ("max-unlocalized-supplementary-alignment-mapping-quality",
     po::value<int>()->default_value(5),
     "Filter reads with supplementary alignments mapped to unlocalized contigs with mapping quality greater than this")

    ("no-reads-with-tag",
     po::value<std::vector<SamTag>>()->multitoken(),
     "Filter reads with tag, use TAG=VALUE to filter specific tag values")
    
    ("disable-downsampling",
     po::bool_switch()->default_value(false),
     "Disables downsampling")
    
    ("downsample-above",
     po::value<int>()->default_value(1000),
     "Downsample reads in regions where coverage is over this")
    
    ("downsample-target",
     po::value<int>()->default_value(500),
     "Target coverage for the downsampler")
    
    ("use-same-read-profile-for-all-samples",
     po::bool_switch()->default_value(false),
     "Use the same read profile for all samples, rather than generating one per sample")
    ;
    
    po::options_description variant_discovery("Variant discovery");
    variant_discovery.add_options()
    ("variant-discovery-mode",
     po::value<CandidateVariantDiscoveryProtocol>()->default_value(CandidateVariantDiscoveryProtocol::illumina),
     "Protocol to use for candidate variant discovery [ILLUMINA, PACBIO]")
     
    ("disable-denovo-variant-discovery",
     po::bool_switch()->default_value(false),
     "Disable all candidate variant discovery from input reads")
    
    ("disable-pileup-candidate-generator",
     po::bool_switch()->default_value(false),
     "Disable candidate generation from raw read alignments (CIGAR strings)")
    
    ("disable-repeat-candidate-generator",
     po::bool_switch()->default_value(false),
     "Disable candidate generation from adjusted read alignments (CIGAR strings) around tandem repeats")
    
    ("disable-assembly-candidate-generator",
     po::bool_switch()->default_value(false),
     "Enable candidate generation using local re-assembly")
    
    ("source-candidates,c",
     po::value<std::vector<fs::path>>()->multitoken(),
     "Variant file paths containing known variants. These variants will automatically become candidates")
    
    ("source-candidates-file",
     po::value<std::vector<fs::path>>()->multitoken(),
     "Files containing lists of source candidate variant files")
     
    ("min-source-candidate-quality",
     po::value<Phred<double>>(),
     "Only variants with quality above this value are considered for candidate generation")
    
    ("use-filtered-source-candidates",
     po::bool_switch()->default_value(false),
     "Use variants from source VCF records that have been filtered")
    
    ("min-pileup-base-quality",
     po::value<int>()->default_value(20),
     "Only bases with quality above this value are considered for candidate generation")
    
    ("min-supporting-reads",
     po::value<int>()->implicit_value(2),
     "Minimum number of reads that must support a variant if it is to be considered a candidate."
     " By default octopus will automatically determine this value")
     
    ("force-pileup-candidates",
     po::bool_switch()->default_value(false),
     "Include pileup candidate variants discovered from reads that are considered likely to be misaligned")
    
    ("max-variant-size",
     po::value<int>()->default_value(2000),
     "Maximum candidate variant size to consider (in region space)")
    
    ("kmer-sizes",
     po::value<std::vector<int>>()->multitoken()
     ->default_value(std::vector<int> {10, 15, 20}, "10 15 20")->composing(),
     "Kmer sizes to use for local assembly")
    
    ("max-fallback-kmers",
     po::value<int>()->default_value(10),
     "How many local assembly fallback kmer sizes to try if the default sizes fail")
    
    ("fallback-kmer-gap",
     po::value<int>()->default_value(10),
     "Gap size used to generate local assembly fallback kmers")
    
    ("max-assembly-region-size",
     po::value<int>()->default_value(600),
     "The maximum region size that can be used for local assembly")
    
    ("max-assembly-region-overlap",
     po::value<int>()->default_value(200),
     "Maximum number of bases allowed to overlap assembly regions")
    
    ("assemble-all",
     po::bool_switch()->default_value(false),
     "Forces all regions to be assembled")
    
    ("assembler-mask-base-quality",
     po::value<int>()->default_value(10),
     "Aligned bases with quality less than this will be converted to reference before"
     " being inserted into the De Bruijn graph")
    
    ("allow-cycles",
     po::bool_switch()->default_value(false),
     "Allow cyclic assembly graphs")
    
    ("min-kmer-prune",
     po::value<int>()->default_value(2),
     "Minimum number of read observations to keep a kmer in the assembly graph before bubble extraction")
    
    ("max-bubbles",
     po::value<int>()->default_value(30),
     "Maximum number of bubbles to extract from the assembly graph")
    
    ("min-bubble-score",
     po::value<double>()->default_value(2.0),
     "Minimum bubble score that will be extracted from the assembly graph")

    ("allow-strand-biased-candidates",
     po::bool_switch()->default_value(false),
     "Do not account for strand bias when evaluating whether to include a candidate variant")
    
    ("min-candidate-credible-vaf-probability",
     po::value<float>()->default_value(0.5),
     "Minimum probability that pileup candidate variant has frequency above '--min-credible-somatic-frequency'")
    ;
    
    po::options_description haplotype_generation("Haplotype generation");
    haplotype_generation.add_options()
    ("max-haplotypes,x",
     po::value<int>()->default_value(200),
     "Maximum number of candidate haplotypes the caller may consider. If a region contains"
     " more candidate haplotypes than this then filtering is applied")
    
    ("haplotype-holdout-threshold",
     po::value<int>()->default_value(2500),
     "Forces the haplotype generator to temporarily hold out some alleles if the number"
     " of haplotypes in a region exceeds this threshold")
    
    ("haplotype-overflow",
     po::value<int>()->default_value(200000),
     "Regions with more haplotypes than this will be skipped")
    
    ("max-holdout-depth",
     po::value<int>()->default_value(20),
     "Maximum number of holdout attempts the haplotype generator can make before the region"
     " is skipped")
    
    ("extension-level",
     po::value<ExtensionLevel>()->default_value(ExtensionLevel::moderate),
     "Level of haplotype extension [MINIMAL, CONSERVATIVE, MODERATE, AGGRESSIVE, UNLIMITED]")
     
    ("lagging-level",
     po::value<LaggingLevel>()->default_value(LaggingLevel::moderate),
     "Level of haplotype lagging [NONE, CONSERVATIVE, MODERATE, OPTIMISTIC, AGGRESSIVE]")
    
    ("backtrack-level",
     po::value<BacktrackLevel>()->default_value(BacktrackLevel::none),
     "Level of backtracking [NONE, MODERATE, AGGRESSIVE]")
    
    ("min-protected-haplotype-posterior",
     po::value<double>()->default_value(1e-10, "1e-10"),
     "Haplotypes with posterior probability less than this may be pruned from the haplotype tree")
    
    ("dont-protect-reference-haplotype",
     po::bool_switch()->default_value(false),
     "Do not protect the reference haplotype from filtering")
    ;
    
    po::options_description general_variant_calling("Variant calling (general)");
    general_variant_calling.add_options()
    ("caller,C",
     po::value<std::string>()->default_value("population"),
     "Which of the octopus calling models to use")
    
    ("organism-ploidy,P",
     po::value<int>()->default_value(2),
     "All contigs with unspecified ploidies are assumed the organism ploidy")
    
    ("contig-ploidies,p",
     po::value<std::vector<ContigPloidy>>()->multitoken()
     ->default_value(std::vector<ContigPloidy> {
        {boost::none, "Y", 1}, {boost::none, "chrY", 1},
        {boost::none, "MT", 1}, {boost::none, "chrM", 1}}, "Y=1 chrY=1 MT=1 chrM=1")
     ->composing(),
     "Space-separated list of contig (contig=ploidy) or sample contig"
     " (sample:contig=ploidy) ploidies")
    
    ("contig-ploidies-file",
     po::value<fs::path>(),
     "File containing a list of contig (contig=ploidy) or sample contig"
     " (sample:contig=ploidy) ploidies, one per line")
    
    ("min-variant-posterior",
     po::value<Phred<double>>()->default_value(Phred<double> {0.1}),
     "Report variant alleles with posterior probability (phred scale) greater than this")
    
    ("refcall",
     po::bool_switch()->default_value(false),
     "Caller will report reference confidence calls for non-variant positions")
     
    ("refcall-block-merge-quality",
     po::value<Phred<double>>()->default_value(Phred<double> {10.0}),
     "Threshold to merge adjacent refcall positions, set to 0 for positional records")
    
    ("max-refcall-posterior",
     po::value<Phred<double>>(),
     "Maximum allowed posterior probability (QUAL) for reference calls")
    
    ("snp-heterozygosity,z",
     po::value<float>()->default_value(0.001, "0.001"),
     "Germline SNP heterozygosity for the given samples")
    
    ("snp-heterozygosity-stdev",
     po::value<float>()->default_value(0.01, "0.01"),
     "Standard deviation of the germline SNP heterozygosity used for the given samples")
    
    ("indel-heterozygosity,y",
     po::value<float>()->default_value(0.0001, "0.0001"),
     "Germline indel heterozygosity for the given samples")
    
    ("use-uniform-genotype-priors",
    po::bool_switch()->default_value(false),
    "Use a uniform prior model when calculating genotype posteriors")
    
    ("max-genotypes",
     po::value<int>(),
     "Maximum number of genotypes that must be evaluated")
    
    ("max-genotype-combinations",
     po::value<int>(),
     "Maximum number of genotype combinations that can be considered when computing joint"
     " genotype posterior probabilities")
    
    ("use-independent-genotype-priors",
     po::bool_switch()->default_value(false),
     "Use independent genotype priors for joint calling")
    
    ("model-posterior",
     po::value<ModelPosteriorPolicy>()->default_value(ModelPosteriorPolicy::special),
     "Policy for calculating model posteriors for variant calls [ALL, OFF, SPECIAL]")
    
    ("disable-inactive-flank-scoring",
     po::bool_switch()->default_value(false),
     "Disable additional calculation to adjust alignment score to account for inactive candidate variants")
    
    ("dont-model-mapping-quality",
     po::bool_switch()->default_value(false),
     "Don't use read mapping quality information in the haplotype likelihood calculation")
    
    ("sequence-error-model",
     po::value<std::string>()->default_value("PCR-free.HiSeq-2500"),
     "Sequencing error model to use by the haplotyoe likelihood model")
    
    ("max-vb-seeds",
     po::value<int>()->default_value(12),
     "Maximum number of seeds to use for Variational Bayes algorithms")
     
    ("max-indel-errors",
     po::value<int>()->default_value(16),
     "Maximum number of indel errors that must be tolerated for haplotype likelihood calculation")
    
    ("use-wide-hmm-scores",
     po::bool_switch()->default_value(false),
     "Use 32-bits rather than 16-bits for HMM scores")

    ("read-linkage",
     po::value<ReadLinkage>()->default_value(ReadLinkage::paired),
     "Read linkage information to use for calling [NONE, PAIRED, LINKED]")
     
    ("min-phase-score",
     po::value<Phred<double>>()->default_value(Phred<double> {5.0}),
     "Minimum phase score (phred scale) required to report sites as phased")
    
    ("phasing-policy",
     po::value<PhasingPolicy>()->default_value(PhasingPolicy::automatic),
     "Policy for applying phasing algorithm [AUTO, CONSERVATIVE, AGGRESSIVE]")

    ("bad-region-tolerance",
     po::value<BadRegionTolerance>()->default_value(BadRegionTolerance::normal),
     "Tolerance for skipping regions that are considered unlikely to be callable [LOW, NORMAL, HIGH, UNLIMITED]")
    ;
    
    po::options_description cancer("Cancer calling model");
    cancer.add_options()
    ("normal-samples,N",
     po::value<std::vector<std::string>>()->multitoken(),
     "Normal samples - all other samples are considered tumour")
    
    ("max-somatic-haplotypes",
     po::value<int>()->default_value(2),
     "Maximum number of somatic haplotypes that may be considered")
    
    ("somatic-snv-prior",
     po::value<float>()->default_value(1e-04, "0.0001"),
     "Prior probability for an SNV somatic mutation at a given base for this sample")
    
    ("somatic-indel-prior",
     po::value<float>()->default_value(1e-06, "0.000001"),
     "Prior probability for an INDEL somatic mutation at a given position for this sample")
    
    ("min-expected-somatic-frequency",
     po::value<float>()->default_value(0.01, "0.01"),
     "Minimum expected somatic allele frequency in the sample")
    
    ("min-credible-somatic-frequency",
     po::value<float>()->default_value(0.005, "0.005"),
     "Minimum credible somatic allele frequency that will be reported")
    
     ("tumour-germline-concentration",
     po::value<float>()->default_value(1.5, "1.5"),
     "Concentration parameter for germline haplotypes in tumour samples")
     
    ("somatic-credible-mass",
     po::value<float>()->default_value(0.9, "0.9"),
     "Mass of the posterior density to use for evaluating somatic allele frequencies")
    
    ("min-somatic-posterior",
     po::value<Phred<double>>()->default_value(Phred<double> {0.5}),
     "Minimum posterior probability (phred scale) to emit a somatic mutation call")
    
    ("normal-contamination-risk",
     po::value<NormalContaminationRisk>()->default_value(NormalContaminationRisk::low),
     "Risk that the normal sample is contaminated by the tumour [LOW, HIGH]")
    
    ("somatics-only",
     po::bool_switch()->default_value(false),
     "Only emit SOMATIC mutations")
    ;
    
    po::options_description trio("Trio calling model");
    trio.add_options()
    ("maternal-sample,M",
     po::value<std::string>(),
     "Maternal sample")
    
    ("paternal-sample,F",
     po::value<std::string>(),
     "Paternal sample")
    
    ("denovo-snv-prior",
     po::value<float>()->default_value(1.3e-8, "1.3e-8"),
     "Prior probability for an SNV de novo mutation at a given base in the offspring")
    
    ("denovo-indel-prior",
     po::value<float>()->default_value(1e-9, "1e-9"),
     "Prior probability for an INDEL de novo mutation at a given position in the offspring")
    
    ("min-denovo-posterior",
     po::value<Phred<double>>()->default_value(Phred<double> {3}),
     "Minimum posterior probability (phred scale) to emit a de novo mutation call")
    
    ("denovos-only",
     po::bool_switch()->default_value(false),
     "Only emit DENOVO mutations")
    ;
    
    po::options_description polyclone("Polyclone calling model");
    polyclone.add_options()
    ("max-clones",
     po::value<int>()->default_value(5),
     "Maximum number of unique clones to consider")
    
    ("min-clone-frequency",
     po::value<float>()->default_value(0.01, "0.01"),
     "Minimum expected clone frequency in the sample")
    
    ("clone-prior",
     po::value<float>()->default_value(0.1, "0.1"),
     "Prior probability of each clone in the sample")
    
    ("clone-concentration",
     po::value<float>()->default_value(1, "1"),
     "Prior mixture concentration for each clone in the sample")
    ;
    
    po::options_description cell("Cell calling model");
    cell.add_options()
    ("max-copy-loss",
     po::value<int>()->default_value(0),
     "Maximum number of haplotype losses in the phylogeny")
    
    ("max-copy-gain",
     po::value<int>()->default_value(0),
     "Maximum number of haplotype gains in the phylogeny")
    
    ("somatic-cnv-prior",
     po::value<float>()->default_value(1e-5, "1e-5"),
     "Prior probability of a given base in a sample being affected by a CNV")
     
    ("dropout-concentration",
    po::value<float>()->default_value(5, "5"),
    "Allelic dropout concentration parameter (default for all samples)")
    
    ("sample-dropout-concentrations",
     po::value<std::vector<SampleDropoutConcentrationPair>>()->multitoken(),
     "Sample allelic dropout concentration parameter (format SAMPLE=CONCENTRATION")
    
    ("phylogeny-concentration",
     po::value<float>()->default_value(20, "20"),
     "Concentration prior for clones in phylogeny")
    ;
    
    po::options_description call_filtering("Call filtering and annotation");
    call_filtering.add_options()
    ("disable-call-filtering",
     po::bool_switch()->default_value(false),
     "Disable call filtering")
    
    ("filter-expression",
     po::value<std::string>()->default_value("QUAL < 10 | MQ < 10 | MP < 10 | AD < 1 | AF < 0.01 | AFB > 0.25 | SB > 0.98 | BQ < 15 | DP < 1 | ADP < 1"),
     "Boolean expression to use to filter variant calls")
    
    ("somatic-filter-expression",
     po::value<std::string>()->default_value("QUAL < 2 | GQ < 20 | MQ < 30 | SMQ < 40 | SB > 0.9 | SD > 0.9 | BQ < 20 | DP < 3 | ADP < 1 | MF > 0.2 | NC > 1 | FRF > 0.5 | AD < 1 | AF < 0.0001"),
     "Boolean expression to use to filter somatic variant calls")
    
    ("denovo-filter-expression",
     po::value<std::string>()->default_value("QUAL < 50 | PP < 40 | GQ < 20 | MQ < 30 | AD < 1 | AF < 0.1 | AFB > 0.2 | SB > 0.95 | BQ < 20 | DP < 10 | ADP < 1 | DC > 1 | MF > 0.2 | FRF > 0.5 | MP < 30 | MQ0 > 2"),
     "Boolean expression to use to filter somatic variant calls")
    
    ("refcall-filter-expression",
     po::value<std::string>()->default_value("QUAL < 2 | GQ < 20 | MQ < 10 | DP < 10 | MF > 0.2"),
     "Boolean expression to use to filter homozygous reference calls")
    
    ("use-preprocessed-reads-for-filtering",
     po::bool_switch()->default_value(false),
     "Use preprocessed reads, as used for calling, for call filtering")
    
    ("keep-unfiltered-calls",
     po::bool_switch()->default_value(false),
     "Keep a copy of unfiltered calls")
     
    ("annotations",
     po::value<std::vector<std::string>>()->multitoken()->implicit_value(std::vector<std::string> {"active"}, "active")->composing(),
     "Annotations to write to final VCF")
    
     ("aggregate-annotations",
     po::bool_switch()->default_value(false),
     "Aggregate all multi-value annotations into a single value")
    
    ("filter-vcf",
     po::value<fs::path>(),
     "Filter the given Octopus VCF without calling")
    
    ("forest-model",
     po::value<fs::path>(),
     "Trained Ranger random forest model file")
    
    ("somatic-forest-model",
     po::value<fs::path>(),
     "Trained Ranger random forest model file for somatic variants")
    
    ("min-forest-quality",
     po::value<Phred<double>>()->default_value(Phred<double> {3}),
     "Minimum PASSing random forest probability (Phred scale)")

    ("use-germline-forest-for-somatic-normals",
     po::bool_switch()->default_value(false),
     "Use the germline forest model for evaluating somatic variant normal sample genotypes rather than the somatic forest model")
    ;
    
    po::options_description all("Octopus command line options");
    all.add(general).add(read_preprocessing)
    .add(variant_discovery).add(haplotype_generation).add(general_variant_calling)
    .add(cancer).add(trio).add(polyclone).add(cell).add(call_filtering);
    
    po::options_description conditional_help_options("Octopus help command line options");
    // Add call_filtering for optional annotation help.
    // Add general_variant_calling to ensure --refcall is resolved correctly.  
    conditional_help_options.add(general).add(general_variant_calling).add(call_filtering);
    OptionMap vm_init;
    po::store(run(po::command_line_parser(argc, argv).options(conditional_help_options).allow_unregistered()), vm_init);
    
    if (vm_init.count("help") == 1) {
        po::store(run(po::command_line_parser(argc, argv).options(general_variant_calling).allow_unregistered()), vm_init);
        if (vm_init.count("caller") == 1 && !vm_init.at("caller").defaulted()) {
            const auto selected_caller = vm_init.at("caller").as<std::string>();
            validate_caller(vm_init);
            if (selected_caller == "individual") {
                po::options_description individual_options("Octopus command line options (individual)");
                individual_options.add(general).add(read_preprocessing)
                .add(variant_discovery).add(haplotype_generation).add(general_variant_calling).add(call_filtering);
                std::cout << individual_options << std::endl;
            } else if (selected_caller == "trio") {
                po::options_description trio_options("Octopus command line options (trio)");
                trio_options.add(general).add(read_preprocessing)
                .add(variant_discovery).add(haplotype_generation).add(general_variant_calling).add(trio).add(call_filtering);
                std::cout << trio_options << std::endl;
            } else if (selected_caller == "population") {
                po::options_description population_options("Octopus command line options (population)");
                population_options.add(general).add(read_preprocessing)
                .add(variant_discovery).add(haplotype_generation).add(general_variant_calling).add(call_filtering);
                std::cout << population_options << std::endl;
            } else if (selected_caller == "cancer") {
                po::options_description cancer_options("Octopus command line options (cancer)");
                cancer_options.add(general).add(read_preprocessing)
                .add(variant_discovery).add(haplotype_generation).add(general_variant_calling).add(cancer).add(call_filtering);
                std::cout << cancer_options << std::endl;
            } else if (selected_caller == "polyclone") {
                po::options_description polyclone_options("Octopus command line options (polyclone)");
                polyclone_options.add(general).add(read_preprocessing)
                .add(variant_discovery).add(haplotype_generation).add(general_variant_calling).add(polyclone).add(call_filtering);
                std::cout << polyclone_options << std::endl;
            } else if (selected_caller == "cell") {
                po::options_description polyclone_options("Octopus command line options (cell)");
                polyclone_options.add(general).add(read_preprocessing)
                                 .add(variant_discovery).add(haplotype_generation).add(general_variant_calling).add(cell)
                                 .add(call_filtering);
                std::cout << polyclone_options << std::endl;
            } else {
                std::cout << all << std::endl;
            }
        } else {
            std::cout << all << std::endl;
        }
        if (vm_init.count("annotations") == 1) {
            std::cout << "Available annotations:\n" << std::endl;
            csr::print_all_measures_help(std::cout);
        }
        return vm_init;
    }
    
    if (vm_init.count("version") == 1) {
        std::cout << "octopus version " << config::Version << '\n'
                  << "Target: " << config::System.system_processor << ' ' << config::System.system_name << " " << config::System.system_version << '\n'
                  << "SIMD extension: " << config::System.simd_extension << '\n'
                  << "Compiler: " << config::System.compiler_name << ' ' << config::System.compiler_version << '\n'
                  << "Boost: " << config::System.boost_version
                  << std::endl;
        return vm_init;
    }
    
    OptionMap vm;
    
    if (vm_init.count("config") == 1) {
        auto config_path = resolve_path(vm_init.at("config").as<fs::path>(), vm_init);
        parse_config_file(config_path, vm, all);
    }
    
    vm_init.clear();
    po::store(run(po::command_line_parser(argc, argv).options(all)), vm);
    validate(vm);
    po::notify(vm);
    
    return vm;
}

class InvalidWorkingDirectory : public UserError
{
    std::string do_where() const override
    {
        return "get_working_directory";
    }
    
    std::string do_why() const override
    {
        std::ostringstream ss {};
        ss << "The working directory you specified ";
        ss << path_;
        ss << " does not exist";
        return ss.str();
    }
    
    std::string do_help() const override
    {
        return "enter a valid working directory";
    }
    
    fs::path path_;
public:
    InvalidWorkingDirectory(fs::path p) : path_ {std::move(p)} {}
};

fs::path get_working_directory(const OptionMap& options)
{
    if (options.count("working-directory") == 1) {
        auto result = expand_user_path(options.at("working-directory").as<fs::path>());
        if (!fs::exists(result) && !fs::is_directory(result)) {
            throw InvalidWorkingDirectory {result};
        }
        return result;
    }
    return fs::current_path();
}

fs::path resolve_path(const fs::path& path, const OptionMap& options)
{
    return ::octopus::resolve_path(path, get_working_directory(options));
}

namespace {

std::string prepend_dashes(std::string option)
{
    option.insert(0, "--");
    return option;
}

std::string implode(std::vector<std::string> options)
{
    std::transform(std::cbegin(options), std::cend(options), std::begin(options), prepend_dashes);
    static const std::string delim {" | "};
    return utils::join(options, delim);
}

}

class CommandLineError : public UserError
{
public:
    CommandLineError() = default;
    
    CommandLineError(std::string&& why) : why_ {std::move(why)} {}
    
protected:
    std::string why_;

private:
    virtual std::string do_where() const override
    {
        return "parse_options";
    }
    
    virtual std::string do_why() const override
    {
        return why_;
    }
    
    virtual std::string do_help() const override
    {
        return "use the --help command to view required and allowable options";
    }
};

class BadConfigFile : public  CommandLineError
{
public:
    BadConfigFile(fs::path p)
    {
        std::ostringstream ss {};
        ss << "The config file path (" << p << ") given in the option '--config' does not exist";
        why_ = ss.str();
    }
};

class UnknownCommandLineOption : public CommandLineError
{
public:
    UnknownCommandLineOption(std::string option)
    : CommandLineError { "The option you specified '--" + option + "' is not recognised"}
    {}
};

class MissingRequiredCommandLineArguement : public CommandLineError
{
public:
    MissingRequiredCommandLineArguement(std::string option)
    : CommandLineError {"The command line option '--" + option + "' is required but is missing"}
    {}
    
    MissingRequiredCommandLineArguement(std::vector<std::string> options, bool strict = false)
    {
        std::ostringstream ss {};
        if (strict) {
            ss << "One ";
        } else {
            ss << "At least one ";
        }
        ss << "of the command line options '" + implode(options) + "' is required but none are present";
        why_ = ss.str();
    }
};

class InvalidCommandLineOptionValue : public CommandLineError
{
public:
    template <typename T>
    InvalidCommandLineOptionValue(std::string option, T value, std::string reason)
    : CommandLineError {
    "The arguement '" + std::to_string(value) + "' given to option '--" + option
    + "' was rejected as it " + reason
    } {}
};

class ConflictingCommandLineOptions : public CommandLineError
{
public:
    ConflictingCommandLineOptions(std::vector<std::string> conflicts)
    {
        std::ostringstream ss {};
        ss << "the options";
        for (const auto& option : conflicts) {
            ss << " " << option;
        }
        ss << " are mutually exclusive";
        why_ = ss.str();
    }
};

class MissingDependentCommandLineOption : public CommandLineError
{
public:
    MissingDependentCommandLineOption(std::string given, std::string dependent)
    {
        std::ostringstream ss {};
        ss << "The option " << given << " requires option " << dependent;
        why_ = ss.str();
    }
};

void parse_config_file(const fs::path& config_file, OptionMap& vm, const po::options_description& options)
{
    if (!fs::exists(config_file)) {
        throw BadConfigFile {config_file};
    }
    std::ifstream config {config_file.string()};
    if (config) {
        try {
            po::store(po::parse_config_file(config, options), vm);
        } catch (const po::invalid_config_file_syntax& e) {
            throw CommandLineError {e.what()};
        } catch (const po::unknown_option& e) {
            throw UnknownCommandLineOption {po::strip_prefixes(e.get_option_name())};
        } catch (const po::invalid_option_value& e) {
            throw CommandLineError {e.what()};
        } catch (const po::invalid_bool_value& e) {
            throw CommandLineError {e.what()};
        } catch (const po::ambiguous_option& e) {
            throw CommandLineError {e.what()};
        } catch (const po::reading_file& e) {
            throw CommandLineError {e.what()};
        }
    }
}

void check_positive(const std::string& option, const OptionMap& vm)
{
    if (vm.count(option) == 1) {
        const auto value = vm.at(option).as<int>();
        if (value < 0) {
            throw InvalidCommandLineOptionValue {option, value, "must be positive" };
        }
    }
}

void check_strictly_positive(const std::string& option, const OptionMap& vm)
{
    if (vm.count(option) == 1) {
        const auto value = vm.at(option).as<int>();
        if (value < 1) {
            throw InvalidCommandLineOptionValue {option, value, "must be greater than zero" };
        }
    }
}

void check_probability(const std::string& option, const OptionMap& vm)
{
    if (vm.count(option) == 1) {
        const auto value = vm.at(option).as<float>();
        if (value < 0 || value > 1) {
            throw InvalidCommandLineOptionValue {option, value, "must be between zero and one" };
        }
    }
}

void conflicting_options(const OptionMap& vm, const std::string& opt1, const std::string& opt2)
{
    if (vm.count(opt1) == 1 && !vm[opt1].defaulted() && vm.count(opt2) == 1 && !vm[opt2].defaulted()) {
        throw ConflictingCommandLineOptions {{opt1, opt2}};
    }
}

void option_dependency(const OptionMap& vm, const std::string& given, const std::string& dependent)
{
    if (vm.count(given) == 1 && !vm[given].defaulted())
        if (vm.count(dependent) == 0 || vm[dependent].defaulted()) {
            throw MissingDependentCommandLineOption {given, dependent};
        }
}

void check_reads_present(const OptionMap& vm)
{
    if (vm.count("reads") == 0 && vm.count("reads-file") == 0) {
        throw MissingRequiredCommandLineArguement {std::vector<std::string> {"reads", "reads-file"}};
    }
}

void check_region_files_consistent(const OptionMap& vm)
{
    if (vm.count("regions-file") == 1 && vm.count("skip-regions-file") == 1) {
        const auto regions_file = vm.at("regions-file").as<fs::path>();
        const auto skip_regions_file = vm.at("skip-regions-file").as<fs::path>();
        if (regions_file == skip_regions_file) {
            throw std::invalid_argument {"options 'regions-file' and 'skip-regions-file' must be unique"};
        }
    }
}

void check_trio_consistent(const OptionMap& vm)
{
    if (vm.at("caller").as<std::string>() == "trio"
        && (vm.count("maternal-sample") == 0 || vm.count("paternal-sample") == 0)) {
        throw std::logic_error {"option 'maternal-sample' and 'paternal-sample' are required"
            " when caller=trio"};
    }
}

void validate_caller(const OptionMap& vm)
{
    if (vm.count("caller") == 1) {
        const auto caller = vm.at("caller").as<std::string>();
        static const std::array<std::string, 6> validCallers {
            "individual", "population", "cancer", "trio", "polyclone", "cell"
        };
        if (std::find(std::cbegin(validCallers), std::cend(validCallers), caller) == std::cend(validCallers)) {
            throw po::validation_error {po::validation_error::kind_t::invalid_option_value, caller, "caller"};
        }
    }
}

po::parsed_options run(po::command_line_parser& parser)
{
    try {
        return parser.run();
    } catch (const po::required_option& e) {
        throw MissingRequiredCommandLineArguement {po::strip_prefixes(e.get_option_name())};
    } catch (const po::unknown_option& e) {
        throw UnknownCommandLineOption {po::strip_prefixes(e.get_option_name())};
    } catch (const po::invalid_option_value& e) {
        throw CommandLineError {e.what()};
    } catch (const po::invalid_bool_value& e) {
        throw CommandLineError {e.what()};
    } catch (const po::ambiguous_option& e) {
        throw CommandLineError {e.what()};
    } catch (const po::reading_file& e) {
        throw CommandLineError {e.what()};
    } catch (const po::invalid_command_line_syntax& e) {
        throw CommandLineError {e.what()};
    } catch (const po::error& e) {
        throw CommandLineError {e.what()};
    }
}

void validate(const OptionMap& vm)
{
    const std::vector<std::string> positive_int_options {
        "threads", "mask-low-quality-tails", "mask-tails", "soft-clip-mask-threshold", "mask-soft-clipped-boundary-bases",
        "min-mapping-quality", "good-base-quality", "min-good-bases", "min-read-length",
        "max-read-length", "min-base-quality", "max-variant-size",
        "max-fallback-kmers", "max-assembly-region-overlap", "assembler-mask-base-quality",
        "min-kmer-prune", "max-bubbles", "max-holdout-depth", "max-copy-loss", "max-copy-gain"
    };
    const std::vector<std::string> strictly_positive_int_options {
        "max-open-read-files", "downsample-above", "downsample-target", "min-supporting-reads",
        "max-assembly-region-size", "fallback-kmer-gap", "organism-ploidy",
        "max-haplotypes", "haplotype-holdout-threshold", "haplotype-overflow",
        "max-genotypes", "max-genotype-combinations", "max-somatic-haplotypes", "max-clones",
        "max-vb-seeds", "max-indel-errors", "max-base-quality", "max-phylogeny-size"
    };
    const std::vector<std::string> probability_options {
        "snp-heterozygosity", "snp-heterozygosity-stdev", "indel-heterozygosity",
        "somatic-snv-prior", "somatic-indel-prior", "min-expected-somatic-frequency", "min-credible-somatic-frequency", "somatic-credible-mass",
        "denovo-snv-prior", "denovo-indel-prior", "min-candidate-credible-vaf-probability",
        "somatic-cnv-prior", "clone-prior"
    };
    conflicting_options(vm, "maternal-sample", "normal-sample");
    conflicting_options(vm, "paternal-sample", "normal-sample");
    conflicting_options(vm, "sites-only", "bamout");
    conflicting_options(vm, "sites-only", "data-profile");
    for (const auto& option : positive_int_options) {
        check_positive(option, vm);
    }
    for (const auto& option : strictly_positive_int_options) {
        check_strictly_positive(option, vm);
    }
    for (const auto& option : probability_options) {
        check_probability(option, vm);
    }
    check_reads_present(vm);
    check_region_files_consistent(vm);
    check_trio_consistent(vm);
    validate_caller(vm);
}

std::istream& operator>>(std::istream& in, ContigPloidy& plodies)
{
    static const std::regex re {"(?:([^:]*):)?([^=]+)=(\\d+)"};
    
    std::string token;
    in >> token;
    std::smatch match;
    
    if (std::regex_match(token, match, re) && match.size() == 4) {
        if (match.length(1) > 0) {
            plodies.sample = match.str(1);
        }
        plodies.contig = match.str(2);
        plodies.ploidy = boost::lexical_cast<decltype(plodies.ploidy)>(match.str(3));
    } else {
        using Error = po::validation_error;
        throw Error {Error::kind_t::invalid_option_value, token, "contig-ploidies"};
    }
    
    return in;
}

std::ostream& operator<<(std::ostream& out, const ContigPloidy& plodies)
{
    if (plodies.sample) out << *plodies.sample << ':';
    out << plodies.contig << "=" << plodies.ploidy;
    return out;
}

std::istream& operator>>(std::istream& in, ContigOutputOrder& order)
{
    std::string token;
    in >> token;
    if (token == "LEXICOGRAPHICAL_ASCENDING")
        order = ContigOutputOrder::lexicographicalAscending;
    else if (token == "LEXICOGRAPHICAL_DESCENDING")
        order = ContigOutputOrder::lexicographicalDescending;
    else if (token == "CONTIG_SIZE_ASCENDING")
        order = ContigOutputOrder::contigSizeAscending;
    else if (token == "CONTIG_SIZE_DESCENDING")
        order = ContigOutputOrder::contigSizeDescending;
    else if (token == "REFERENCE_INDEX")
        order = ContigOutputOrder::referenceIndex;
    else if (token == "REFERENCE_INDEX_REVERSED")
        order = ContigOutputOrder::referenceIndexReversed;
    else if (token == "UNSPECIFIED")
        order = ContigOutputOrder::unspecified;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "contig-output-order"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ContigOutputOrder& order)
{
    switch (order) {
        case ContigOutputOrder::lexicographicalAscending:
            out << "LEXICOGRAPHICAL_ASCENDING";
            break;
        case ContigOutputOrder::lexicographicalDescending:
            out << "LEXICOGRAPHICAL_DESCENDING";
            break;
        case ContigOutputOrder::contigSizeAscending:
            out << "CONTIG_SIZE_ASCENDING";
            break;
        case ContigOutputOrder::contigSizeDescending:
            out << "CONTIG_SIZE_DESCENDING";
            break;
        case ContigOutputOrder::referenceIndex:
            out << "REFERENCE_INDEX";
            break;
        case ContigOutputOrder::referenceIndexReversed:
            out << "REFERENCE_INDEX_REVERSED";
            break;
        case ContigOutputOrder::unspecified:
            out << "UNSPECIFIED";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ExtensionLevel& level)
{
    std::string token;
    in >> token;
    if (token == "MINIMAL")
        level = ExtensionLevel::minimal;
    else if (token == "CONSERVATIVE")
        level = ExtensionLevel::conservative;
    else if (token == "MODERATE")
        level = ExtensionLevel::moderate;
    else if (token == "AGGRESSIVE")
        level = ExtensionLevel::aggressive;
    else if (token == "UNLIMITED")
        level = ExtensionLevel::unlimited;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "extension-level"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ExtensionLevel& level)
{
    switch (level) {
    case ExtensionLevel::minimal:
        out << "MINIMAL";
        break;
    case ExtensionLevel::conservative:
        out << "CONSERVATIVE";
        break;
    case ExtensionLevel::moderate:
        out << "MODERATE";
        break;
    case ExtensionLevel::aggressive:
        out << "AGGRESSIVE";
        break;
    case ExtensionLevel::unlimited:
        out << "UNLIMITED";
        break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, BacktrackLevel& level)
{
    std::string token;
    in >> token;
    if (token == "NONE")
        level = BacktrackLevel::none;
    else if (token == "MODERATE")
        level = BacktrackLevel::moderate;
    else if (token == "AGGRESSIVE")
        level = BacktrackLevel::aggressive;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "backtrack-level"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const BacktrackLevel& level)
{
    switch (level) {
    case BacktrackLevel::none:
        out << "NONE";
        break;
    case BacktrackLevel::moderate:
        out << "MODERATE";
        break;
    case BacktrackLevel::aggressive:
        out << "AGGRESSIVE";
        break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, LaggingLevel& level)
{
    std::string token;
    in >> token;
    if (token == "NONE")
        level = LaggingLevel::none;
    else if (token == "CONSERVATIVE")
        level = LaggingLevel::conservative;
    else if (token == "MODERATE")
        level = LaggingLevel::moderate;
    else if (token == "OPTIMISTIC")
        level = LaggingLevel::optimistic;
    else if (token == "AGGRESSIVE")
        level = LaggingLevel::aggressive;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "lagging-level"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const LaggingLevel& level)
{
    switch (level) {
        case LaggingLevel::none:
            out << "NONE";
            break;
        case LaggingLevel::conservative:
            out << "CONSERVATIVE";
            break;
        case LaggingLevel::moderate:
            out << "MODERATE";
            break;
        case LaggingLevel::optimistic:
            out << "OPTIMISTIC";
            break;
        case LaggingLevel::aggressive:
            out << "AGGRESSIVE";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, NormalContaminationRisk& result)
{
    std::string token;
    in >> token;
    if (token == "LOW")
        result = NormalContaminationRisk::low;
    else if (token == "HIGH")
        result = NormalContaminationRisk::high;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "normal-contamination-risk"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const NormalContaminationRisk& risk)
{
    switch (risk) {
        case NormalContaminationRisk::low:
            out << "LOW";
            break;
        case NormalContaminationRisk::high:
            out << "HIGH";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, BadRegionTolerance& result)
{
    std::string token;
    in >> token;
    if (token == "LOW")
        result = BadRegionTolerance::low;
    else if (token == "NORMAL")
        result = BadRegionTolerance::normal;
    else if (token == "HIGH")
        result = BadRegionTolerance::high;
    else if (token == "UNLIMITED")
        result = BadRegionTolerance::unlimited;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "bad-region-tolerance"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const BadRegionTolerance& tolerance)
{
    switch (tolerance) {
        case BadRegionTolerance::low:
            out << "LOW";
            break;
        case BadRegionTolerance::normal:
            out << "NORMAL";
            break;
        case BadRegionTolerance::high:
            out << "HIGH";
            break;
        case BadRegionTolerance::unlimited:
            out << "UNLIMITED";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ReadLinkage& result)
{
    std::string token;
    in >> token;
    if (token == "NONE")
        result = ReadLinkage::none;
    else if (token == "PAIRED")
        result = ReadLinkage::paired;
    else if (token == "LINKED")
        result = ReadLinkage::linked;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "read-linkage"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ReadLinkage& linkage)
{
    switch (linkage) {
        case ReadLinkage::none:
            out << "NONE";
            break;
        case ReadLinkage::paired:
            out << "PAIRED";
            break;
        case ReadLinkage::linked:
            out << "LINKED";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, CandidateVariantDiscoveryProtocol& result)
{
    std::string token;
    in >> token;
    if (token == "ILLUMINA")
        result = CandidateVariantDiscoveryProtocol::illumina;
    else if (token == "PACBIO")
        result = CandidateVariantDiscoveryProtocol::pacbio;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "variant-discovery-protocol"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const CandidateVariantDiscoveryProtocol& protocol)
{
    switch (protocol) {
        case CandidateVariantDiscoveryProtocol::illumina:
            out << "ILLUMINA";
            break;
        case CandidateVariantDiscoveryProtocol::pacbio:
            out << "PACBIO";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, RealignedBAMType& result)
{
    std::string token;
    in >> token;
    if (token == "FULL")
        result = RealignedBAMType::full;
    else if (token == "MINI")
        result = RealignedBAMType::mini;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "bamout-type"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const RealignedBAMType& type)
{
    switch (type) {
        case RealignedBAMType::full:
            out << "FULL";
            break;
        case RealignedBAMType::mini:
            out << "MINI";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ReadDeduplicationDetectionPolicy& result)
{
    std::string token;
    in >> token;
    if (token == "RELAXED")
        result = ReadDeduplicationDetectionPolicy::relaxed;
    else if (token == "AGGRESSIVE")
        result = ReadDeduplicationDetectionPolicy::aggressive;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "bamout-type"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ReadDeduplicationDetectionPolicy& policy)
{
    switch (policy) {
        case ReadDeduplicationDetectionPolicy::relaxed:
            out << "RELAXED";
            break;
        case ReadDeduplicationDetectionPolicy::aggressive:
            out << "AGGRESSIVE";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ModelPosteriorPolicy& result)
{
    std::string token;
    in >> token;
    if (token == "ALL")
        result = ModelPosteriorPolicy::all;
    else if (token == "OFF")
        result = ModelPosteriorPolicy::off;
    else if (token == "SPECIAL")
        result = ModelPosteriorPolicy::special;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "model-posterior"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ModelPosteriorPolicy& policy)
{
    switch (policy) {
        case ModelPosteriorPolicy::all:
            out << "ALL";
            break;
        case ModelPosteriorPolicy::off:
            out << "OFF";
            break;
        case ModelPosteriorPolicy::special:
            out << "SPECIAL";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, SampleDropoutConcentrationPair& result)
{
    std::string token;
    in >> token;
    const auto equal_pos = token.find_last_of("=");
    const auto concentration = token.substr(equal_pos + 1);
    token.erase(equal_pos);
    result.sample = std::move(token);
    result.concentration = boost::lexical_cast<float>(concentration);
    return in;
}

std::ostream& operator<<(std::ostream& os, const SampleDropoutConcentrationPair& concentration)
{
    os << concentration.sample << '=' << concentration.concentration;
    return os;
}

std::istream& operator>>(std::istream& in, PhasingPolicy& result)
{
    std::string token;
    in >> token;
    if (token == "CONSERVATIVE")
        result = PhasingPolicy::conservative;
    else if (token == "AGGRESSIVE")
        result = PhasingPolicy::aggressive;
    else if (token == "AUTO")
        result = PhasingPolicy::automatic;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token, "phasing-policy"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const PhasingPolicy& policy)
{
    switch (policy) {
        case PhasingPolicy::conservative:
            out << "CONSERVATIVE";
            break;
        case PhasingPolicy::aggressive:
            out << "AGGRESSIVE";
            break;
        case PhasingPolicy::automatic:
            out << "AUTO";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, SamTag& result)
{
    std::string token;
    in >> token;
    const auto eq_n = token.find_last_of('=');
    if (eq_n != std::string::npos) {
        result.value = token.substr(eq_n + 1);
        token.erase(eq_n);
    }
    result.tag = token;
    return in;
}

std::ostream& operator<<(std::ostream& out, const SamTag& tag)
{
    out << tag.tag;
    if (tag.value) out << '=' << *tag.value;
    return out;
}

namespace {

template <typename T>
bool is_type(const OptionMap::mapped_type& value)
{
    try {
        boost::any_cast<T>(value.value());
        return true;
    } catch (const boost::bad_any_cast&) {
        return false;
    }
}

template <typename T>
bool is_vector_type(const OptionMap::mapped_type& value)
{
    return is_type<std::vector<T>>(value);
}

bool is_list_type(const OptionMap::mapped_type& value)
{
    return is_vector_type<int>(value) || is_vector_type<float>(value) || is_vector_type<double>(value)
            || is_vector_type<std::string>(value) || is_vector_type<fs::path>(value) || is_vector_type<ContigPloidy>(value);
}

template <typename T>
void write_vector(const OptionMap& options, const OptionMap::key_type& label, std::ostream& os, const char bullet)
{
    const auto& values = options[label].as<std::vector<T>>();
    std::size_t i {0};
    for (const auto& val : values) {
        os << bullet << ' ' << label << "[" << i++ << "]=" << val;
        if (i != values.size()) os << std::endl;
    }
}

template <>
void write_vector<fs::path>(const OptionMap& options, const OptionMap::key_type& label, std::ostream& os, const char bullet)
{
    const auto& values = options[label].as<std::vector<fs::path>>();
    std::size_t i {0};
    for (const auto& val : values) {
        os << bullet << ' ' << label << "[" << i++ << "]=" << val.filename().string();
        if (i != values.size()) os << std::endl;
    }
}

} // namespace

std::ostream& operator<<(std::ostream& os, const OptionMap& options)
{
    std::size_t i {0};
    for (const auto& p : options) {
        const auto& label = p.first;
        const auto& value = p.second;
        const char bullet {options[label].defaulted() || value.defaulted() ? '>' : '~'};
        if (!is_list_type(value)) {
            os << bullet << ' ' << label;
            if (((boost::any)value.value()).empty()) {
                os << "(empty)";
            }
            os << "=";
        }
        if (static_cast<boost::any>(value.value()).type() == typeid(int)) {
            os << options[label].as<int>();
        } else if (static_cast<boost::any>(value.value()).type() == typeid(bool)) {
            os << (options[label].as<bool>() ? "yes" : "no");
        } else if (static_cast<boost::any>(value.value()).type() == typeid(float)) {
            os << options[label].as<float>();
        } else if (static_cast<boost::any>(value.value()).type() == typeid(double)) {
            os << options[label].as<double>();
        } else if (is_type<const char*>(value)) {
            os << options[label].as<const char*>();
        } else if (is_type<std::string>(value)) {
            const auto str = options[label].as<std::string>();
            if (str.size()) {
                os << str;
            } else {
                os << "true";
            }
        } else if (is_type<fs::path>(value)) {
            os << options[label].as<fs::path>().filename().string();
        } else if (is_vector_type<int>(value)) {
            write_vector<int>(options, label, os, bullet);
        } else if (is_vector_type<float>(value)) {
            write_vector<float>(options, label, os, bullet);
        } else if (is_vector_type<double>(value)) {
            write_vector<double>(options, label, os, bullet);
        } else if (is_vector_type<std::string>(value)) {
            write_vector<std::string>(options, label, os, bullet);
        } else if (is_vector_type<fs::path>(value)) {
            write_vector<fs::path>(options, label, os, bullet);
        } else if (is_type<Phred<double>>(value)) {
            os << options[label].as<Phred<double>>();
        } else if (is_type<MemoryFootprint>(value)) {
            os << options[label].as<MemoryFootprint>();
        } else if (is_type<ContigOutputOrder>(value)) {
            os << options[label].as<ContigOutputOrder>();
        } else if (is_type<ContigPloidy>(value)) {
            os << options[label].as<ContigPloidy>();
        } else if (is_vector_type<ContigPloidy>(value)) {
            write_vector<ContigPloidy>(options, label, os, bullet);
        } else if (is_type<ExtensionLevel>(value)) {
            os << options[label].as<ExtensionLevel>();
        } else if (is_type<LaggingLevel>(value)) {
            os << options[label].as<LaggingLevel>();
        } else if (is_type<BacktrackLevel>(value)) {
            os << options[label].as<BacktrackLevel>();
        } else if (is_type<NormalContaminationRisk>(value)) {
            os << options[label].as<NormalContaminationRisk>();
        } else if (is_type<BadRegionTolerance>(value)) {
            os << options[label].as<BadRegionTolerance>();
        } else if (is_type<ReadLinkage>(value)) {
            os << options[label].as<ReadLinkage>();
        } else if (is_type<CandidateVariantDiscoveryProtocol>(value)) {
            os << options[label].as<CandidateVariantDiscoveryProtocol>();
        } else if (is_type<RealignedBAMType>(value)) {
            os << options[label].as<RealignedBAMType>();
        } else if (is_type<ReadDeduplicationDetectionPolicy>(value)) {
            os << options[label].as<ReadDeduplicationDetectionPolicy>();
        } else if (is_vector_type<SampleDropoutConcentrationPair>(value)) {
            write_vector<SampleDropoutConcentrationPair>(options, label, os, bullet);
        } else if (is_type<ModelPosteriorPolicy>(value)) {
            os << options[label].as<ModelPosteriorPolicy>();
        } else if (is_type<PhasingPolicy>(value)) {
            os << options[label].as<PhasingPolicy>();
        } else if (is_type<SamTag>(value)) {
            os << options[label].as<SamTag>();
        } else if (is_vector_type<SamTag>(value)) {
            write_vector<SamTag>(options, label, os, bullet);
        } else {
            os << "UnknownType(" << ((boost::any)value.value()).type().name() << ")";
        }
        if (++i != options.size()) os << std::endl;
    }
    return os;
}

std::string to_string(const OptionMap& options, const bool one_line, const bool mark_modified)
{
    std::ostringstream ss {};
    ss << options;
    auto result = ss.str();
    if (one_line) {
        auto chunks = utils::split(result, '\n');
        for (auto& chunk : chunks) {
            const bool modified {chunk[0] == '~'};
            chunk[0] = '-';
            chunk[1] = '-';
            if (modified && mark_modified) {
                chunk.insert(0, 1,'*');
            }
            chunk[chunk.find('=')] = ' ';
        }
        result = utils::join(chunks, ' ');
    }
    return result;
}

} // namespace options
} // namespace octopus
