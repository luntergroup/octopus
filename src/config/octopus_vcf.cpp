// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "octopus_vcf.hpp"

#include <unordered_map>
#include <functional>

namespace octopus { namespace vcf {

VcfHeader::Builder make_header_template()
{
    VcfHeader::Builder result {};
    
    result.set_file_format(vcfspec::version);
    
    result.add_info("AA", "1", "String", "Ancestral allele");
    result.add_info("AC", "A", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed");
    //result.add_info("AF", "A", "Float", "Allele Frequency, for each ALT allele, in the same order as listed");
    result.add_info("AN", "1", "Integer", "Total number of alleles in called genotypes");
    result.add_info("BQ", "1", "Integer", "RMS base quality at this position");
    result.add_info("DP", "1", "Integer", "Combined depth across samples");
    //result.add_info("END", "1", "Integer", "End position of the variant described in this record");
    result.add_info("MQ", "1", "Integer", "RMS mapping quality");
    result.add_info("MQ0", "1", "Integer", "Number of MAPQ == 0 reads covering this record");
    result.add_info("NS", "1", "Integer", "Number of samples with data");
    result.add_info("SB", "1", "Float", "Strand bias at this position");
    
    result.add_format("GT", "1", "String", "Genotype");
    result.add_format("DP", "1", "Integer", "Read depth at this position for this sample");
    //result.add_format("FT", "1", "String", "Sample genotype filter indicating if this genotype was “called”");
    //result.add_format("GL", "G", "Float", "log10-scaled genotype likelihoods");
    //result.add_format("GLE", "1", "Integer", "Genotype likelihoods of heterogeneous ploidy");
    //result.add_format("PL", "G", "Integer", "Phred-scaled genotype likelihoods");
    //result.add_format("GP", "G", "Float", "Phred-scaled genotype posterior probabilities");
    result.add_format("GQ", "1", "Integer", "Conditional genotype quality (phred-scaled)");
    result.add_format("PS", "1", "String", "Phase set");
    result.add_format("PQ", "1", "Integer", "Phasing quality");
    result.add_format("MQ", "1", "Integer", "RMS mapping quality");
    result.add_format("BQ", "1", "Integer", "RMS base quality at this position");
    
    result.add_filter("PASS", "All filters passed");
    
    return result;
}

static const std::unordered_map<std::string, std::string> filter_descriptions
{
{spec::filter::q10, "Variant quality is below 10"},
{spec::filter::q20, "Variant quality is below 20"},
{spec::filter::lowQuality, "Variant quality is low"},
{spec::filter::highMappingQualityDivergence, "High Kullback–Leibler divergence between REF and ALT mapping quality distributions"},
{spec::filter::alleleBias, "The called allele frequencies are not as expected for the given ploidy"},
{spec::filter::lowModelPosterior, "Variant failed model posterior filter"},
{spec::filter::lowMappingQuality, "Mapping quality across calling region is low"},
{spec::filter::highMappingQualityZeroCount, "The number of reads with mapping quality zero is low"},
{spec::filter::lowQualityByDepth, "Variant failed quality/depth filter"},
{spec::filter::strandBias, "One of the alternative alleles has strand bias"},
{spec::filter::lowDepth, "Read depth around variant is low"},
{spec::filter::filteredReadFraction, "The number of reads filtered for calling is high"},
{spec::filter::lowGQ, "Sample genotype quality low"}
};

VcfHeader::Builder& add_filter(VcfHeader::Builder& builder, const std::string& key)
{
    builder.add_filter(key, filter_descriptions.at(key));
    return builder;
}

} // namespace vcf
} // namespace octopus
