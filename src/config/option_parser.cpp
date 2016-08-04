// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "option_parser.hpp"

#include <regex>
#include <iostream>
#include <stdexcept>

#include <basics/phred.hpp>

#include "common.hpp"

namespace po = boost::program_options;

namespace octopus { namespace options {

void conflicting_options(const OptionMap& vm, const std::string& opt1, const std::string& opt2);
void option_dependency(const OptionMap& vm, const std::string& for_what, const std::string& required_option);
void check_reads_present(const OptionMap& vm);
void check_region_files_consistent(const OptionMap& vm);
void check_trio_consistent(const OptionMap& vm);
void validate_caller(const OptionMap& vm);
void validate_options(const OptionMap& vm);

OptionMap parse_options(const int argc, const char** argv)
{
    po::options_description general("General");
    general.add_options()
    ("help,h", "Produce help message")
    
    ("version", "Output the version number")
    
    ("config", "A config file, used to populate command line options")
    
    ("debug",
     po::bool_switch()->default_value(false),
     "Writes verbose debug information to debug.log in the working directory")
    
    ("trace",
     po::bool_switch()->default_value(false),
     "Writes very verbose debug information to trace.log in the working directory")
    ;
    
    po::options_description backend("Backend");
    backend.add_options()
    ("working-directory,wd",
     po::value<std::string>(),
     "Sets the working directory")
    
    ("threads,t",
     po::value<unsigned>()->implicit_value(0),
     "Maximum number of threads to be used, enabling this option with no argument lets the application"
     " decide the number of threads ands enables specific algorithm parallelisation")
    
    ("max-reference-cache-footprint,mrcf",
     po::value<float>()->default_value(50),
     "Maximum memory footprint for cached reference sequence (in megabytes)")
    
    ("target-read-buffer-footprint,trbf",
     po::value<float>()->default_value(0.5),
     "None binding request to limit the memory footprint of buffered read data (in gigabytes)")
    
    ("compress-reads,cr",
     po::bool_switch()->default_value(false),
     "Compresses all read data when not being used resulting in a smaller memory footprint"
     " but slower processing")
    
    ("max-open-read-files,morf",
     po::value<unsigned>()->default_value(250),
     "Limits the number of read files that can be open simultaneously")
    ;
    
    po::options_description input("I/O");
    input.add_options()
    ("reference,R",
     po::value<std::string>()->required(),
     "FASTA format reference genome file to be analysed. Target regions"
     " will be extracted from the reference index if not provded explicitly")
    
    ("reads,r",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-seperated list of BAM/CRAM files to be analysed."
     " May be specified multiple times")
    
    ("reads-file,rf",
     po::value<std::string>(),
     "File containing a list of BAM/CRAM files, one per line, to be analysed")
    
    ("one-based-indexing,1bi",
     po::bool_switch()->default_value(false),
     "Notifies that input regions are given using one based indexing rather than zero based")
    
    ("regions,T",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-seperated list of regions (chrom:begin-end) to be analysed."
     " May be specified multiple times")
    
    ("regions-file,TF",
     po::value<std::string>(),
     "File containing a list of regions (chrom:begin-end), one per line, to be analysed")
    
    ("skip-regions,sr",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-seperated list of regions (chrom:begin-end) to skip"
     " May be specified multiple times")
    
    ("skip-regions-file,srf",
     po::value<std::string>(),
     "File of regions (chrom:begin-end), one per line, to skip")
    
    ("samples,S",
     po::value<std::vector<std::string>>()->multitoken(),
     "Space-seperated list of sample names to analyse")
    
    ("samples-file,SF",
     po::value<std::string>(),
     "File of sample names to analyse, one per line, which must be a subset of the samples"
     " that appear in the read files")
    
    ("output,o",
     po::value<std::string>()->default_value("octopus_calls.vcf"),
     "File to where output is written")
    
    ("contig-output-order,coo",
     po::value<ContigOutputOrder>()->default_value(ContigOutputOrder::AsInReferenceIndex),
     "The order contigs should be written to the output")
    
    ("legacy",
     po::bool_switch()->default_value(false),
     "Outputs a legacy version of the final callset in addition to the native version")
    
    ("regenotype",
     po::value<std::string>(),
     "VCF file specifying calls to regenotype, only sites in this files will appear in the"
     " final output")
    ;
    
    po::options_description transforms("Read transformations");
    transforms.add_options()
    ("disable-read-transforms",
     po::bool_switch()->default_value(false),
     "Disables all read transformations")
    
    ("disable-soft-clip-masking",
     po::bool_switch()->default_value(false),
     "Disables soft clipped masking, thus allowing all soft clipped bases to be used"
     " for candidate generation")
    
    ("mask-tails",
     po::value<unsigned>()->implicit_value(3),
     "Masks this number of bases of the tail of all reads")
    
    ("mask-soft-clipped-boundries",
     po::value<unsigned>()->default_value(2),
     "Masks this number of adjacent non soft clipped bases when soft clipped bases are present")
    
    ("disable-adapter-masking",
     po::bool_switch()->default_value(false),
     "Disables adapter detection and masking")
    
    ("disable-overlap-masking",
     po::bool_switch()->default_value(false),
     "Disables read segment overlap masking")
    ;
    
    po::options_description filters("Read filtering");
    filters.add_options()
    ("disable-read-filtering",
     po::bool_switch()->default_value(false),
     "Disables all read filters")
    
    ("consider-unmapped-reads,allow-unmapped",
     po::bool_switch()->default_value(false),
     "Allows reads marked as unmapped to be used for calling")
    
    ("min-mapping-quality,min-mq",
     po::value<unsigned>()->default_value(20),
     "Minimum read mapping quality required to consider a read for calling")
    
    ("good-base-quality,good-bq",
     po::value<unsigned>()->default_value(20),
     "Base quality threshold used by min-good-bases and min-good-base-fraction filters")
    
    ("min-good-base-fraction,min-good-bp-frac",
     po::value<double>()->implicit_value(0.5),
     "Base quality threshold used by min-good-bases filter")
    
    ("min-good-bases,min-good-bps",
     po::value<unsigned>()->default_value(20),
     "Minimum number of bases with quality min-base-quality before read is considered")
    
    ("allow-qc-fails",
     po::bool_switch()->default_value(false),
     "Filters reads marked as QC failed")
    
    ("min-read-length,min-read-len",
     po::value<unsigned>(),
     "Filters reads shorter than this")
    
    ("max-read-length,max-read-len",
     po::value<unsigned>(),
     "Filter reads longer than this")
    
    ("allow-marked-duplicates,allow-marked-dups",
     po::bool_switch()->default_value(false),
     "Allows reads marked as duplicate in alignment record")
    
    ("allow-octopus-duplicates,allow-dups",
     po::bool_switch()->default_value(false),
     "Allows reads considered duplicates by octopus")
    
    ("no-secondary-alignments",
     po::bool_switch()->default_value(false),
     "Filters reads marked as secondary alignments")
    
    ("no-supplementary-alignmenets",
     po::bool_switch()->default_value(false),
     "Filters reads marked as supplementary alignments")
    
    ("consider-reads-with-unmapped-segments",
     po::bool_switch()->default_value(false),
     "Allows reads with unmapped template segmenets to be used for calling")
    
    ("consider-reads-with-distant-segments",
     po::bool_switch()->default_value(false),
     "Allows reads with template segmenets that are on different contigs")
    
    ("allow-adapter-contaminated-reads",
     po::bool_switch()->default_value(false),
     "Allows reads with possible adapter contamination")
    
    ("disable-downsampling,no-downsampling",
     po::bool_switch()->default_value(false),
     "Diables all downsampling")
    
    ("downsample-above",
     po::value<unsigned>()->default_value(500),
     "Downsample reads in regions where coverage is over this")
    
    ("downsample-target",
     po::value<unsigned>()->default_value(400),
     "The target coverage for the downsampler")
    ;
    
    po::options_description candidates("Candidate variant generation");
    candidates.add_options()
    ("disable-raw-cigar-candidate-generator,no-cigar-candidates",
     po::bool_switch()->default_value(false),
     "Disables candidate generation from raw read alignments (CIGAR strings)")
    
    ("disable-assembly-candidate-generator,no-assembly-candidates",
     po::bool_switch()->default_value(false),
     "Disables candidate generation using local re-assembly")
    
    ("generate-candidates-from-source,source",
     po::value<std::string>(),
     "Variant file path containing known variants. These variants will automatically become"
     " candidates")
    
    ("min-base-quality,min-bq",
     po::value<unsigned>()->default_value(20),
     "Only bases with quality above this value are considered for candidate generation")
    
    ("min-supporting-reads,min-support",
     po::value<unsigned>()->implicit_value(2),
     "Minimum number of reads that must support a variant if it is to be considered a candidate."
     " By default octopus will automatically determine this value")
    
    ("max-variant-size,max-var-size",
     po::value<unsigned>()->default_value(2000),
     "Maximum candidate varaint size to consider (in region space)")
    
    ("kmer-size,kmer",
     po::value<std::vector<unsigned>>()->multitoken()
     ->default_value(std::vector<unsigned> {10, 25}, "10 25")->composing(),
     "K-mer sizes to use for local re-assembly")
    
    ("assembler-mask-base-quality",
     po::value<unsigned>()->implicit_value(10),
     "Matching alignment bases with quality less than this will be reference masked before."
     " Ff no value is specified then min-base-quality is used")
    ;
    
    po::options_description caller("Common caller options");
    caller.add_options()
    ("caller,C",
     po::value<std::string>()->default_value("population"),
     "Which of the octopus callers to use")
    
    ("organism-ploidy,ploidy",
     po::value<unsigned>()->default_value(2),
     "All contigs with unspecified ploidies are assumed the organism ploidy")
    
    ("contig-ploidies",
     po::value<std::vector<ContigPloidy>>()->multitoken()
     ->default_value(std::vector<ContigPloidy> {
        {boost::none, "Y", 1}, {boost::none, "MT", 1}}, "Y=1 MT=1")
     ->composing(),
     "Space-seperated list of contig (contig=ploidy) or sample contig"
     " (sample:contig=ploidy) ploidies")
    
    ("contig-ploidies-file",
     po::value<std::string>(),
     "File containing a list of contig (contig=ploidy) or sample contig"
     " (sample:contig=ploidy) ploidies, one per line")
    
    ("min-variant-posterior,min-post",
     po::value<Phred<double>>()->default_value(Phred<double> {2.0}),
     "Report variant alleles with posterior probability (phred scale) greater than this")
    
    ("min-refcall-posterior,min-ref-post",
     po::value<Phred<double>>()->default_value(Phred<double> {2.0}),
     "Report reference alleles with posterior probability (phred scale) greater than this")
    
    ("report-refcalls,gvcf",
     po::value<RefCallType>()->implicit_value(RefCallType::Blocked),
     "Caller will report reference confidence calls for each position (Positional),"
     " or in automatically sized blocks (Blocked)")
    
    ("sites-only",
     po::bool_switch()->default_value(false),
     "Only reports call sites (i.e. without sample genotype information)")
    
    ("snp-heterozygosity,snp-hets",
     po::value<float>()->default_value(0.001, "0.001"),
     "The germline SNP heterozygosity used to calculate genotype priors")
    
    ("indel-heterozygosity,indel-hets",
     po::value<float>()->default_value(0.0001, "0.0001"),
     "The germline indel heterozygosity used to calculate genotype priors")
    ;
    
    po::options_description cancer("Cancer caller");
    cancer.add_options()
    ("normal-sample,normal",
     po::value<std::string>(),
     "Normal sample - all other samples are considered tumour")
    
    ("somatic-mutation-rate,somatic-rate",
     po::value<float>()->default_value(0.00001, "0.00001"),
     "Expected somatic mutation rate, per megabase pair, for this sample")
    
    ("min-somatic-frequency,min-somatic-freq",
     po::value<float>()->default_value(0.01, "0.01"),
     "minimum allele frequency that can be considered as a viable somatic mutation")
    
    ("credible-mass,cm",
     po::value<float>()->default_value(0.99, "0.99"),
     "Mass of the posterior density to use for evaluating allele frequencies")
    
    ("min-somatic-posterior,min-somatic-post",
     po::value<Phred<double>>()->default_value(Phred<double> {2.0}),
     "Minimum somatic mutation call posterior probability (phred scale)")
    
    ("somatics-only",
     po::bool_switch()->default_value(false),
     "Only report somatic variant calls")
    ;
    
    po::options_description trio("Trio caller");
    trio.add_options()
    ("maternal-sample,mother",
     po::value<std::string>(),
     "Maternal sample")
    
    ("paternal-sample,father",
     po::value<std::string>(),
     "Paternal sample")
    
    ("denovos-only",
     po::bool_switch()->default_value(false),
     "Only report de novo variant calls (i.e. alleles unique to the child)")
    ;
    
    po::options_description phaser("Phasing options");
    phaser.add_options()
    ("phasing-level,phase",
     po::value<PhasingLevel>()->default_value(PhasingLevel::Conservative),
     "Level of phasing - longer range phasing can improve calling accuracy at the cost"
     " of runtime speed. Possible values are: Minimal, Conservative, Aggressive")
    
    ("min-phase-score",
     po::value<Phred<double>>()->default_value(Phred<double> {20.0}),
     "Minimum phase score (phred scale) required to report sites as phased")
    
    ("use-unconditional-phase-score",
     po::bool_switch()->default_value(false),
     "Computes unconditional phase scores rather than conditioning on called genotypes")
    
    ("disable-read-guided-phasing",
     po::bool_switch()->default_value(false),
     "Restricts phase score computation to use only genotype posteriors")
    ;
    
    po::options_description advanced("Advanced calling algorithm");
    advanced.add_options()
    ("max-haplotypes,max-haps",
     po::value<unsigned>()->default_value(128),
     "Maximum number of candidate haplotypes the caller may consider")
    
    ("min-haplotype-filter-posterior,min-hap-post",
     po::value<float>()->default_value(1e-10, "1e-10"),
     "Haplotypes with posterior less than this can be filtered, allowing greater"
     " longer haplotype extesion in complex regions")
    
    ("disable-inactive-flank-scoring,noIFS",
     po::bool_switch()->default_value(false),
     "Disables additional calculation to adjust alignment score when there are inactive"
     " candidates in haplotype flanking regions")
    ;
    
    po::options_description call_filtering("Callset filtering");
    call_filtering.add_options()
    ("disable-call-filtering,no-filtering",
     po::bool_switch()->default_value(false),
     "Disables all callset filtering")
    
    ("disable-model-filtering,noMF",
     po::bool_switch()->default_value(false),
     "Disables model based filtering of variant calls")
    ;
    
    po::options_description all("octopus options");
    all.add(general).add(backend).add(input).add(transforms).add(filters)
    .add(candidates).add(caller).add(advanced).add(cancer).add(trio).add(phaser)
    .add(call_filtering);
    
    OptionMap vm_init;
    
    po::store(po::command_line_parser(argc, argv).options(general)
              .allow_unregistered().run(), vm_init);
    
    if (vm_init.count("help") == 1) {
        po::store(po::command_line_parser(argc, argv).options(caller)
                  .allow_unregistered().run(), vm_init);
        
        if (vm_init.count("caller") == 1) {
            const auto caller = vm_init.at("caller").as<std::string>();
            
            validate_caller(vm_init);
            
            if (caller == "individual") {
                std::cout << all << std::endl;
            } else if (caller == "population") {
                std::cout << all << std::endl;
            } else if (caller == "cancer") {
                std::cout << all << std::endl;
            } else {
                std::cout << all << std::endl;
            }
        } else {
            std::cout << all << std::endl;
        }
        
        return vm_init;
    }
    
    if (vm_init.count("version") == 1) {
        std::cout << "octopus " << info::VERSION << std::endl;
        return vm_init;
    }
    
    OptionMap vm;
    
    if (vm_init.count("config") == 1) {
        std::ifstream config {vm.at("config").as<std::string>()};
        if (config) {
            po::store(po::parse_config_file(config, all), vm);
        }
    }
    
    vm_init.clear();
    
    po::store(po::command_line_parser(argc, argv).options(all).run(), vm);
    
    // boost::option cannot handle option dependencies so we must do our own checks
    validate_options(vm);
    
    po::notify(vm);
    
    return vm;
}

void conflicting_options(const OptionMap& vm, const std::string& opt1, const std::string& opt2)
{
    if (vm.count(opt1) && !vm[opt1].defaulted() && vm.count(opt2) && !vm[opt2].defaulted()) {
        throw std::logic_error(std::string("conflicting options '") + opt1 + "' and '" + opt2 + "'.");
    }
}

void option_dependency(const OptionMap& vm, const std::string& for_what,
                       const std::string& required_option)
{
    if (vm.count(for_what) && !vm[for_what].defaulted())
        if (vm.count(required_option) == 0 || vm[required_option].defaulted()) {
            throw std::logic_error(std::string("option '") + for_what
                                   + "' requires option '" + required_option + "'.");
        }
}

void check_reads_present(const OptionMap& vm)
{
    if (vm.count("reads") == 0 && vm.count("reads-file") == 0) {
        throw po::required_option {"--reads | --reads-file"};
    }
}

void check_region_files_consistent(const OptionMap& vm)
{
    if (vm.count("regions-file") == 1 && vm.count("skip-regions-file") == 1) {
        const auto regions_file = vm.at("regions-file").as<std::string>();
        const auto skip_regions_file = vm.at("skip-regions-file").as<std::string>();
        if (regions_file == skip_regions_file) {
            throw std::invalid_argument {"options 'regions-file' and 'skip-regions-file' must"
                " have unique values"};
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
        
        static const std::vector<std::string> valid_callers {
            "individual", "population", "cancer", "trio"
        };
        
        if (std::find(std::cbegin(valid_callers), std::cend(valid_callers), caller)
            == std::cend(valid_callers)) {
            throw po::validation_error {po::validation_error::kind_t::invalid_option_value, caller,
                "caller"};
        }
    }
}

void validate_options(const OptionMap& vm)
{
    check_reads_present(vm);
    check_region_files_consistent(vm);
    check_trio_consistent(vm);
    validate_caller(vm);
}
    
std::istream& operator>>(std::istream& in, ContigPloidy& result)
{
    static const std::regex re {"(?:([^:]*):)?([^=]+)=(\\d+)"};
    
    std::string token;
    in >> token;
    
    std::smatch match;
    if (std::regex_match(token, match, re) && match.size() == 4) {
        if (match.length(1) > 0) {
            result.sample = match.str(1);
        }
        
        result.contig = match.str(2);
        result.ploidy = boost::lexical_cast<decltype(result.ploidy)>(match.str(3));
    } else {
        using Error = po::validation_error;
        throw Error {Error::kind_t::invalid_option_value, token, "contig-ploidies"};
    }
    
    return in;
}

std::ostream& operator<<(std::ostream& out, const ContigPloidy& cp)
{
    if (cp.sample) out << *cp.sample << ':';
    out << cp.contig << "=" << cp.ploidy;
    return out;
}

std::istream& operator>>(std::istream& in, RefCallType& result)
{
    std::string token;
    in >> token;
    if (token == "Positional")
        result = RefCallType::Positional;
    else if (token == "Blocked")
        result = RefCallType::Blocked;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token,
        "refcalls"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const RefCallType& type)
{
    switch (type) {
        case RefCallType::Positional:
            out << "Positional";
            break;
        case RefCallType::Blocked:
            out << "Blocked";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, ContigOutputOrder& result)
{
    std::string token;
    in >> token;
    if (token == "LexicographicalAscending")
        result = ContigOutputOrder::LexicographicalAscending;
    else if (token == "LexicographicalDescending")
        result = ContigOutputOrder::LexicographicalDescending;
    else if (token == "ContigSizeAscending")
        result = ContigOutputOrder::ContigSizeAscending;
    else if (token == "ContigSizeDescending")
        result = ContigOutputOrder::ContigSizeDescending;
    else if (token == "AsInReference")
        result = ContigOutputOrder::AsInReferenceIndex;
    else if (token == "AsInReferenceReversed")
        result = ContigOutputOrder::AsInReferenceIndexReversed;
    else if (token == "Unspecified")
        result = ContigOutputOrder::Unspecified;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token,
        "contig-output-order"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const ContigOutputOrder& order)
{
    switch (order) {
        case ContigOutputOrder::LexicographicalAscending:
            out << "LexicographicalAscending";
            break;
        case ContigOutputOrder::LexicographicalDescending:
            out << "LexicographicalDescending";
            break;
        case ContigOutputOrder::ContigSizeAscending:
            out << "ContigSizeAscending";
            break;
        case ContigOutputOrder::ContigSizeDescending:
            out << "ContigSizeDescending";
            break;
        case ContigOutputOrder::AsInReferenceIndex:
            out << "AsInReferenceIndex";
            break;
        case ContigOutputOrder::AsInReferenceIndexReversed:
            out << "AsInReferenceIndexReversed";
            break;
        case ContigOutputOrder::Unspecified:
            out << "Unspecified";
            break;
    }
    return out;
}

std::istream& operator>>(std::istream& in, PhasingLevel& result)
{
    std::string token;
    in >> token;
    if (token == "Minimal")
        result = PhasingLevel::Minimal;
    else if (token == "Conservative")
        result = PhasingLevel::Conservative;
    else if (token == "Aggressive")
        result = PhasingLevel::Aggressive;
    else throw po::validation_error {po::validation_error::kind_t::invalid_option_value, token,
        "phasing-level"};
    return in;
}

std::ostream& operator<<(std::ostream& out, const PhasingLevel& level)
{
    switch (level) {
        case PhasingLevel::Minimal:
            out << "Minimal";
            break;
        case PhasingLevel::Conservative:
            out << "Conservative";
            break;
        case PhasingLevel::Aggressive:
            out << "Aggressive";
            break;
    }
    return out;
}
} // namespace options
} // namespace octopus
