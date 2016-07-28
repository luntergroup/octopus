//
//  quality_by_depth.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "quality_by_depth.hpp"

#include "vcf_record.hpp"

namespace octopus { namespace CallFiltering
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
} // namespace CallFiltering
} // namespace octopus
