// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef octopus_vcf_hpp
#define octopus_vcf_hpp

#include <string>

#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_spec.hpp"

namespace octopus { namespace vcf {

namespace spec {

namespace info {

VCF_SPEC_CONSTANT modelPosterior {"MP"};
VCF_SPEC_CONSTANT denovo {"DENOVO"};
VCF_SPEC_CONSTANT reversion {"REVERSION"};

} // namespace info

namespace filter {

VCF_SPEC_CONSTANT q3 {"q3"};
VCF_SPEC_CONSTANT q5 {"q5"};
VCF_SPEC_CONSTANT q10 {"q10"};
VCF_SPEC_CONSTANT q20 {"q20"};
VCF_SPEC_CONSTANT lowQuality {"LQ"};
VCF_SPEC_CONSTANT lowPosteriorProbability {"LPP"};
VCF_SPEC_CONSTANT lowDepth {"DP"};
VCF_SPEC_CONSTANT highMappingQualityDivergence {"MQD"};
VCF_SPEC_CONSTANT alleleFrequencyBias {"AFB"};
VCF_SPEC_CONSTANT lowAlleleFrequency {"AF"};
VCF_SPEC_CONSTANT alleleDepth {"AD"};
VCF_SPEC_CONSTANT lowModelPosterior {"MP"};
VCF_SPEC_CONSTANT lowMappingQuality {"MQ"};
VCF_SPEC_CONSTANT highMappingQualityZeroCount {"MQ0"};
VCF_SPEC_CONSTANT lowQualityByDepth {"QD"};
VCF_SPEC_CONSTANT strandBias {"SB"};
VCF_SPEC_CONSTANT filteredReadFraction {"FRF"};
VCF_SPEC_CONSTANT highGCRegion {"GC"};
VCF_SPEC_CONSTANT lowGQ {"GQ"};
VCF_SPEC_CONSTANT lowGQD {"GQD"};
VCF_SPEC_CONSTANT highClippedReadFraction {"CRF"};
VCF_SPEC_CONSTANT bq10 {"bq10"};
VCF_SPEC_CONSTANT lowBaseQuality {"LBQ"};
VCF_SPEC_CONSTANT highMismatchCount {"MC"};
VCF_SPEC_CONSTANT highMismatchFraction {"MF"};
VCF_SPEC_CONSTANT normalContamination {"NC"};
VCF_SPEC_CONSTANT deNovoContamination {"DC"};
VCF_SPEC_CONSTANT readSideBias {"RSB"};
VCF_SPEC_CONSTANT readTailBias {"RTB"};
VCF_SPEC_CONSTANT readEndBias {"REB"};
VCF_SPEC_CONSTANT strandDisequilibrium {"SD"};
VCF_SPEC_CONSTANT lowClassificationConfidence {"CC"};
VCF_SPEC_CONSTANT somaticMappingQuality {"SMQ"};
VCF_SPEC_CONSTANT lowAssignedDepth {"ADP"};
VCF_SPEC_CONSTANT lowDuplicateConcordance {"LDC"};
VCF_SPEC_CONSTANT highDuplicateAlleleDepth {"DAD"};
VCF_SPEC_CONSTANT highDuplicateAlleleFraction {"DAF"};
VCF_SPEC_CONSTANT highErrorRate {"ER"};
VCF_SPEC_CONSTANT highErrorRateStdev {"ERS"};

} // namespace filter

namespace format {

VCF_SPEC_CONSTANT modelPosterior {"MP"};

} // namespace info

} // namespace spec

VcfHeader::Builder make_header_template();

VcfHeader::Builder& add_filter(VcfHeader::Builder& builder, const std::string& key);

} // namespace vcf
} // namespace octopus

#endif
