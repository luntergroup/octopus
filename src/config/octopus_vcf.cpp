// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "octopus_vcf.hpp"

#include <unordered_map>
#include <functional>

namespace octopus { namespace vcf {

VcfHeader::Builder make_header_template()
{
    VcfHeader::Builder result {};
    
    result.set_file_format(vcfspec::version);
    
    result.add_info("AC", "A", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed");
    result.add_info("AN", "1", "Integer", "Total number of alleles in called genotypes");
    result.add_info("DP", "1", "Integer", "Combined depth across samples");
    result.add_info("MQ", "1", "Float", "RMS mapping quality");
    result.add_info("MQ0", "1", "Integer", "Number of MAPQ == 0 reads covering this record");
    result.add_info("NS", "1", "Integer", "Number of samples with data");
    result.add_info("END", "1", "Integer", "End position on CHROM");
    
    result.add_format("GT", "1", "String", "Genotype");
    result.add_format("DP", "1", "Integer", "Read depth at this position for this sample");
    result.add_format("FT", "1", "String", "Sample genotype filter indicating if this genotype was 'called'");
    result.add_format("GQ", "1", "Integer", "Conditional genotype quality (phred-scaled)");
    result.add_format("PS", "1", "Integer", "Phase set");
    result.add_format("PQ", "1", "Integer", "Phasing quality");
    result.add_format("MQ", "1", "Integer", "RMS mapping quality");
    
    result.add_filter("PASS", "All filters passed");
    
    return result;
}

static const std::unordered_map<std::string, std::string> filter_descriptions
{
{spec::filter::q3, "Variant quality is below 3"},
{spec::filter::q5, "Variant quality is below 5"},
{spec::filter::q10, "Variant quality is below 10"},
{spec::filter::q20, "Variant quality is below 20"},
{spec::filter::lowQuality, "Variant quality is low"},
{spec::filter::lowPosteriorProbability, "Variant posterior probability is low"},
{spec::filter::highMappingQualityDivergence, "High Kullback–Leibler divergence between REF and ALT mapping quality distributions"},
{spec::filter::lowAlleleFrequency, "The emperical allele frequency is low"},
{spec::filter::alleleFrequencyBias, "The called allele frequencies are not as expected for the given ploidy"},
{spec::filter::lowModelPosterior, "Variant failed model posterior filter"},
{spec::filter::lowMappingQuality, "Mapping quality across calling region is low"},
{spec::filter::highMappingQualityZeroCount, "The number of reads with mapping quality zero is low"},
{spec::filter::lowQualityByDepth, "Variant failed quality/depth filter"},
{spec::filter::strandBias, "One of the alternative alleles has strand bias"},
{spec::filter::lowDepth, "Read depth around variant is low"},
{spec::filter::filteredReadFraction, "The number of reads filtered for calling is high"},
{spec::filter::highGCRegion, "The GC content of the region is too high"},
{spec::filter::lowGQ, "Sample genotype quality low"},
{spec::filter::lowGQD, "Sample GQD is low"},
{spec::filter::highClippedReadFraction, "High fraction of clipped reads covering position"},
{spec::filter::bq10, "Median base quality supporting variant is less than 10"},
{spec::filter::lowBaseQuality, "Median base quality supporting variant is low"},
{spec::filter::highMismatchFraction, "Count of reads containing mismatch to called allele is high"},
{spec::filter::highMismatchFraction, "Fraction of reads containing mismatch to called allele is high"},
{spec::filter::normalContamination, "Somatic haplotypes detected in a called normal sample"},
{spec::filter::deNovoContamination, "De novo allele detected in the offsprings parents"},
{spec::filter::readSideBias, "Variant is biased to one side of supporting reads"},
{spec::filter::readTailBias, "Variant is biased towards tail of supporting reads"},
{spec::filter::readEndBias, "Variant is biased towards ends (head or tail) of supporting reads"},
{spec::filter::strandDisequilibrium, "Reads overlapping the site are biased to one strand"},
{spec::filter::lowClassificationConfidence, "Classification confidence is low"},
{spec::filter::alleleDepth, "Low empirical allele depth"},
{spec::filter::somaticMappingQuality, "Low somatic mapping quality"},
{spec::filter::lowAssignedDepth, "Assigned depth is low"},
{spec::filter::lowDuplicateConcordance, "Duplicate read concordance is low"},
{spec::filter::highDuplicateAlleleDepth, "Number of duplicate reads supporting allele is high"},
{spec::filter::highDuplicateAlleleFraction, "High fraction of reads supporting allele are duplicates"},
{spec::filter::highErrorRate, "High error rate in reads overlapping the site"},
{spec::filter::highErrorRateStdev, "High error rate standard deviation in reads overlapping the site"}
};

VcfHeader::Builder& add_filter(VcfHeader::Builder& builder, const std::string& key)
{
    builder.add_filter(key, filter_descriptions.at(key));
    return builder;
}

} // namespace vcf
} // namespace octopus
