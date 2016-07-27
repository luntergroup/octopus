//
//  quality_by_depth.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef quality_by_depth_hpp
#define quality_by_depth_hpp

#include <string>
#include <memory>
#include <utility>

#include "measure.hpp"

class VcfRecord;

namespace Octopus { namespace CallFiltering
{
    class QualityByDepth : public Measure
    {
        virtual double operator()(const VcfRecord& call) const override;
        virtual std::string name() const override;
    };
} // namespace CallFiltering
} // namespace Octopus

#endif /* quality_by_depth_hpp */
