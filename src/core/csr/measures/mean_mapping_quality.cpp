// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mean_mapping_quality.hpp"

#include "io/variant/vcf_record.hpp"

namespace octopus { namespace csr {

double MeanMappingQuality::operator()(const VcfRecord& call) const
{
    return std::stod(call.info_value("MQ").front());
}

std::string MeanMappingQuality::name() const
{
    return "mmq";
}

} // namespace csr
} // namespace octopus
