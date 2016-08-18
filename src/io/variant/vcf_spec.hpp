// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef vcf_spec_h
#define vcf_spec_h

#include <array>

#define VCF_SPEC_SEPERATOR static constexpr char

#define VCF_SPEC_CONSTANT static constexpr const char*

/**
 The namespace vcfspec contains all the reserved definitions used in the VCF specification.
 See vcfspec::Version for the VCF version the remaining definitions conform to.
 */
namespace octopus { namespace vcfspec {

    VCF_SPEC_CONSTANT version {"VCFv4.3"};
    
    VCF_SPEC_CONSTANT missingValue {"."};
    
    namespace header {
    
    VCF_SPEC_CONSTANT lineOpener {"#"};
    
    VCF_SPEC_CONSTANT chrom {"CHROM"};
    
    VCF_SPEC_CONSTANT pos {"POS"};
    
    VCF_SPEC_CONSTANT id {"ID"};
    
    VCF_SPEC_CONSTANT ref {"REF"};
    
    VCF_SPEC_CONSTANT alt {"ALT"};
    
    VCF_SPEC_CONSTANT qual {"QUAL"};
    
    VCF_SPEC_CONSTANT filter {"FILTER"};
    
    VCF_SPEC_CONSTANT info {"INFO"};
    
    VCF_SPEC_CONSTANT format {"FORMAT"};
    
    static constexpr std::array<const char*, 8> requiredFields {
        chrom, pos, id, ref, alt, qual, filter, info
    };
    
    static constexpr std::array<const char*, 9> requiredFieldsWithSamples {
        chrom, pos, id, ref, alt, qual, filter, info, format
    };
    
    namespace meta {
    
    namespace tag {
    
    VCF_SPEC_CONSTANT info {"INFO"};
    
    VCF_SPEC_CONSTANT filter {"FILTER"};
    
    VCF_SPEC_CONSTANT format {"FORMAT"};
    
    VCF_SPEC_CONSTANT contig {"contig"};
    
    } // namespace tag
    
    namespace struc {
    
    VCF_SPEC_CONSTANT id {"ID"};
    
    VCF_SPEC_CONSTANT number {"Number"};
    
    VCF_SPEC_CONSTANT type {"Type"};
    
    VCF_SPEC_CONSTANT description {"Description"};
    
    VCF_SPEC_CONSTANT source {"Source"};
    
    VCF_SPEC_CONSTANT version {"Version"};
    
    // Not all these fields are required in a structured field, but those present
    // should be in this order
    static constexpr std::array<const char*, 6> order {
        id, number, type, description, source, version
    };
        
    } // namespace struc
    
    VCF_SPEC_CONSTANT vcfVersion {"fileformat"};
    
    } // namespace meta
    
    } // namespace header
    
    namespace allele {
    
    VCF_SPEC_SEPERATOR seperator {','};
    
    } // namespace allele
    
    namespace filter {

    VCF_SPEC_SEPERATOR seperator {';'};
    
    VCF_SPEC_CONSTANT pass {"PASS"};
    
    } // namespace filter
    
    namespace info {
    
    VCF_SPEC_SEPERATOR seperator {';'};
    
    VCF_SPEC_SEPERATOR valueSeperator {','};
    
    VCF_SPEC_CONSTANT ancestralAllele {"AA"};
    
    VCF_SPEC_CONSTANT alleleReadDepth {"AD"};
    
    VCF_SPEC_CONSTANT alleleReadDepthForward {"ADF"};
    
    VCF_SPEC_CONSTANT alleleReadDepthReverse {"ADR"};
        
    VCF_SPEC_CONSTANT combinedReadDepth {"DP"};
    
    VCF_SPEC_CONSTANT alleleFrequency {"AF"};
    
    VCF_SPEC_CONSTANT totalAlleleCount {"AN"};
    
    VCF_SPEC_CONSTANT rmsBaseQuality {"BQ"};
    
    VCF_SPEC_CONSTANT dbSNPMember {"DB"};
    
    VCF_SPEC_CONSTANT endPosition {"END"};
    
    VCF_SPEC_CONSTANT hapmap2Member {"H2"};
    
    VCF_SPEC_CONSTANT hapmap3Member {"H3"};
    
    VCF_SPEC_CONSTANT rmsMappingQuality {"MQ"};
    
    VCF_SPEC_CONSTANT mappingQualityZeroCount {"MQ0"};
    
    VCF_SPEC_CONSTANT hapmap3 {"H3"};
    
    VCF_SPEC_CONSTANT numSamplesWithData {"NS"};
    
    VCF_SPEC_CONSTANT strandBias {"SB"};
    
    VCF_SPEC_CONSTANT somatic {"SOMATIC"};
    
    VCF_SPEC_CONSTANT validated {"VALIDATED"};
    
    VCF_SPEC_CONSTANT thousandGenomes {"1000G"};
        
    } // namespace info
    
    namespace format {
    
    VCF_SPEC_SEPERATOR seperator {':'};
    
    VCF_SPEC_SEPERATOR valueSeperator {','};
    
    VCF_SPEC_CONSTANT genotype {"GT"};
    
    VCF_SPEC_CONSTANT alleleReadDepth {"AD"};
    
    VCF_SPEC_CONSTANT alleleReadDepthForward {"ADF"};
    
    VCF_SPEC_CONSTANT alleleReadDepthReverse {"ADR"};
    
    VCF_SPEC_CONSTANT combinedReadDepth {"DP"};
    
    VCF_SPEC_CONSTANT filter {"FT"};
    
    VCF_SPEC_CONSTANT conditionalQuality {"GQ"};
    
    VCF_SPEC_CONSTANT posteriors {"GP"};
    
    VCF_SPEC_CONSTANT likelihoods {"GL"};
    
    VCF_SPEC_CONSTANT haplotypeQualities {"HQ"};
    
    VCF_SPEC_CONSTANT rmsMappingQuality {"MQ"};
    
    VCF_SPEC_CONSTANT phredLikelihoods {"PL"};
    
    VCF_SPEC_CONSTANT phaseQuality {"PQ"};
    
    VCF_SPEC_CONSTANT phaseSet {"PS"};
        
    VCF_SPEC_CONSTANT phasedSeperator {"|"};
    
    VCF_SPEC_CONSTANT unphasedSeperator {"/"};
    
    } // namespace format
    
} // namespace vcfspec
} // namespace octopus

#endif
