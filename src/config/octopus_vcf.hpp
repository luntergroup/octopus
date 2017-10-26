// Copyright (c) 2016 Daniel Cooke
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

} // namespace info

namespace filter {

VCF_SPEC_CONSTANT q10 {"q10"};
VCF_SPEC_CONSTANT q20 {"q20"};
VCF_SPEC_CONSTANT mappingQualityDivergence {"MQD"};
VCF_SPEC_CONSTANT alleleBias {"AFB"};
VCF_SPEC_CONSTANT modelPosterior {"MP"};
VCF_SPEC_CONSTANT mappingQuality {"MQ"};
VCF_SPEC_CONSTANT qualityByDepth {"QD"};
VCF_SPEC_CONSTANT strandBias {"SB"};

} // namespace filter

} // namespace spec

VcfHeader::Builder make_header_template();

VcfHeader::Builder& add_filter(VcfHeader::Builder& builder, const std::string& key);

} // namespace vcf
} // namespace octopus

#endif
