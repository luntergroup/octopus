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

    VCF_SPEC_CONSTANT Version {"VCFv4.3"};
    
    VCF_SPEC_CONSTANT MissingValue {"."};
    
    namespace header {
    
    VCF_SPEC_CONSTANT LineOpener {"#"};
    
    VCF_SPEC_CONSTANT Chrom {"CHROM"};
    
    VCF_SPEC_CONSTANT Pos {"POS"};
    
    VCF_SPEC_CONSTANT Id {"ID"};
    
    VCF_SPEC_CONSTANT Ref {"REF"};
    
    VCF_SPEC_CONSTANT Alt {"ALT"};
    
    VCF_SPEC_CONSTANT Qual {"QUAL"};
    
    VCF_SPEC_CONSTANT Filter {"FILTER"};
    
    VCF_SPEC_CONSTANT Info {"INFO"};
    
    VCF_SPEC_CONSTANT Format {"FORMAT"};
    
    static constexpr std::array<const char*, 8> RequiredFields {
        Chrom, Pos, Id, Ref, Alt, Qual, Filter, Info
    };
    
    static constexpr std::array<const char*, 9> RequiredFieldsWithSamples {
        Chrom, Pos, Id, Ref, Alt, Qual, Filter, Info, Format
    };
    
    namespace meta {
    
    namespace tag {
    
    VCF_SPEC_CONSTANT Info {"INFO"};
    
    VCF_SPEC_CONSTANT Filter {"FILTER"};
    
    VCF_SPEC_CONSTANT Format {"FORMAT"};
    
    VCF_SPEC_CONSTANT Contig {"contig"};
    
    } // namespace tag
    
    namespace struc {
    
    VCF_SPEC_CONSTANT Id {"ID"};
    
    VCF_SPEC_CONSTANT Number {"Number"};
    
    VCF_SPEC_CONSTANT Type {"Type"};
    
    VCF_SPEC_CONSTANT Description {"Description"};
    
    VCF_SPEC_CONSTANT Source {"Source"};
    
    VCF_SPEC_CONSTANT Version {"Version"};
    
    // Not all these fields are required in a structured field, but those present
    // should be in this order
    static constexpr std::array<const char*, 6> Order {
        Id, Number, Type, Description, Source, Version
    };
        
    } // namespace struc
    
    VCF_SPEC_CONSTANT VcfVersion {"fileformat"};
    
    } // namespace meta
    
    } // namespace header
    
    namespace allele {
    
    VCF_SPEC_SEPERATOR Seperator {','};
    
    } // namespace allele
    
    namespace filter {

    static constexpr char Seperator {';'};
    
    VCF_SPEC_CONSTANT Pass {"PASS"};
    
    } // namespace filter
    
    namespace info {
    
    VCF_SPEC_SEPERATOR Seperator {';'};
    
    VCF_SPEC_SEPERATOR ValueSeperator {','};
    
    VCF_SPEC_CONSTANT AncestralAllele {"AA"};
    
    VCF_SPEC_CONSTANT AlleleReadDepth {"AD"};
    
    VCF_SPEC_CONSTANT AlleleReadDepthForward {"ADF"};
    
    VCF_SPEC_CONSTANT AlleleReadDepthReverse {"ADR"};
        
    VCF_SPEC_CONSTANT CombinedReadDepth {"DP"};
    
    VCF_SPEC_CONSTANT AlleleFrequency {"AF"};
    
    VCF_SPEC_CONSTANT TotalAlleleCount {"AN"};
    
    VCF_SPEC_CONSTANT RMSBaseQuality {"BQ"};
    
    VCF_SPEC_CONSTANT DbSNPMember {"DB"};
    
    VCF_SPEC_CONSTANT EndPosition {"END"};
    
    VCF_SPEC_CONSTANT Hapmap2Member {"H2"};
    
    VCF_SPEC_CONSTANT Hapmap3Member {"H3"};
    
    VCF_SPEC_CONSTANT RMSMappingQuality {"MQ"};
    
    VCF_SPEC_CONSTANT MappingQualityZeroCount {"MQ0"};
    
    VCF_SPEC_CONSTANT Hapmap3 {"H3"};
    
    VCF_SPEC_CONSTANT NumSamplesWithData {"NS"};
    
    VCF_SPEC_CONSTANT StrandBias {"SB"};
    
    VCF_SPEC_CONSTANT Somatic {"SOMATIC"};
    
    VCF_SPEC_CONSTANT Validated {"VALIDATED"};
    
    VCF_SPEC_CONSTANT ThousandGenomes {"1000G"};
        
    } // namespace info
    
    namespace format {
    
    VCF_SPEC_SEPERATOR Seperator {':'};
    
    VCF_SPEC_SEPERATOR ValueSeperator {','};
    
    VCF_SPEC_CONSTANT Genotype {"GT"};
    
    VCF_SPEC_CONSTANT AlleleReadDepth {"AD"};
    
    VCF_SPEC_CONSTANT AlleleReadDepthForward {"ADF"};
    
    VCF_SPEC_CONSTANT AlleleReadDepthReverse {"ADR"};
    
    VCF_SPEC_CONSTANT CombinedReadDepth {"DP"};
    
    VCF_SPEC_CONSTANT Filter {"FT"};
    
    VCF_SPEC_CONSTANT ConditionalQuality {"GQ"};
    
    VCF_SPEC_CONSTANT Posteriors {"GP"};
    
    VCF_SPEC_CONSTANT Likelihoods {"GL"};
    
    VCF_SPEC_CONSTANT HaplotypeQualities {"HQ"};
    
    VCF_SPEC_CONSTANT RMSMappingQuality {"MQ"};
    
    VCF_SPEC_CONSTANT PhredLikelihoods {"PL"};
    
    VCF_SPEC_CONSTANT PhaseQuality {"PQ"};
    
    VCF_SPEC_CONSTANT PhaseSet {"PS"};
        
    VCF_SPEC_CONSTANT PhasedSeperator {"|"};
    
    VCF_SPEC_CONSTANT UnphasedSeperator {"/"};
    
    } // namespace format
    
} // namespace vcfspec
} // namespace octopus

#endif
