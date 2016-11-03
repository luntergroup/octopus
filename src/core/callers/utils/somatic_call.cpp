// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_call.hpp"

#include "utils/string_utils.hpp"

namespace octopus {

void SomaticCall::decorate(VcfRecord::Builder& record) const
{
    record.set_somatic();
    record.add_format("SCR");
    
    for (const auto& p : credible_regions_) {
        if (p.second.somatic) {
            using utils::to_string;
            record.set_format(p.first, "SCR", {
                    to_string(p.second.somatic->first, 2),
                    to_string(p.second.somatic->second, 2)
            });
        } else {
            record.set_format(p.first, "SCR", {"0", "0"});
        }
    }
}

} // namespace octopus
