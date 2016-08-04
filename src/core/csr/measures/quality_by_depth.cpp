// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "quality_by_depth.hpp"

#include <io/variant/vcf_record.hpp>

namespace octopus { namespace csr
{
    double QualityByDepth::operator()(const VcfRecord& call) const
    {
        if (call.qual()) {
            return *call.qual() / std::stod(call.info_value("DP").front());
        } else {
            return 0;
        }
    }
    
    std::string QualityByDepth::name() const
    {
        return "QD";
    }
} // namespace csr
} // namespace octopus
