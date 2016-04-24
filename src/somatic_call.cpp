//
//  somatic_call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "somatic_call.hpp"

namespace Octopus
{
    void SomaticCall::decorate(VcfRecord::Builder& record) const
    {
        record.set_somatic();
    }
} // namespace Octopus
