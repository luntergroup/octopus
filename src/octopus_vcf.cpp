//
//  octopus_vcf.cpp
//  Octopus
//
//  Created by Daniel Cooke on 05/06/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "octopus_vcf.hpp"

namespace octopus
{
    namespace Vcf
    {
        VcfHeader::Builder make_octopus_header_template()
        {
            VcfHeader::Builder result {};
            
            result.add_info("AA", "1", "String", "Ancestral allele");
            result.add_info("AC", "1", "Integer", "Allele count in genotypes, for each ALT allele, in the same order as listed");
            result.add_info("AF", "A", "Float", "Allele Frequency, for each ALT allele, in the same order as listed");
            result.add_info("AN", "1", "Integer", "Total number of alleles in called genotypes");
            result.add_info("BQ", "1", "Integer", "RMS base quality at this position");
            result.add_info("DP", "1", "Integer", "Combined depth across samples");
            result.add_info("END", "1", "Integer", "End position of the variant described in this record");
            result.add_info("MQ", "1", "Integer", "RMS mapping quality");
            result.add_info("MQ0", "1", "Integer", "Number of MAPQ == 0 reads covering this record");
            result.add_info("NS", "1", "Integer", "Number of samples with data");
            result.add_info("SB", "1", "Float", "Strand bias at this position");
            result.add_info("SOMATIC", "0", "Flag", "Indicates that the record is a somatic mutation, for cancer genomics");
            
            result.add_format("GT", "1", "String", "Genotype");
            result.add_format("DP", "1", "Integer", "Read depth at this position for this sample");
            result.add_format("FT", "1", "String", "Sample genotype filter indicating if this genotype was “called”");
            result.add_format("GL", "G", "Float", "log10-scaled genotype likelihoods");
            result.add_format("GLE", "1", "Integer", "Genotype likelihoods of heterogeneous ploidy");
            result.add_format("PL", "G", "Integer", "Phred-scaled genotype likelihoods");
            result.add_format("GP", "G", "Float", "Phred-scaled genotype posterior probabilities");
            result.add_format("GQ", "1", "Integer", "Conditional genotype quality (phred-scaled)");
            result.add_format("PS", "1", "String", "Phase set");
            result.add_format("PQ", "1", "Integer", "Phasing quality");
            result.add_format("MQ", "1", "Integer", "RMS mapping quality");
            result.add_format("BQ", "1", "Integer", "RMS base quality at this position");
            
            result.add_filter("PASS", "All filters passed");
            result.add_filter("MQ", "Root-mean-square mapping quality across calling region is low");
            result.add_filter("q10", "Variant quality is below 10");
            result.add_filter("SB", "One of the alternative alleles has strand bias");
            result.add_filter("KL", "High Kullback–Leibler divergence between REF and ALT mapping quality distributions");
            
            return result;
        }
    } // namespace Vcf
} // namespace octopus
