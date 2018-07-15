// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_call.hpp"

#include "utils/string_utils.hpp"

namespace octopus {

void SomaticCall::decorate(VcfRecord::Builder& record) const
{
    record.set_somatic();
    if (posterior_) {
        record.set_info("PP", utils::to_string(posterior_->score()));
    }
    if (!map_vafs_.empty()) {
        record.add_format("MAP_VAF");
    }
    record.add_format("VAF_CR");
    for (const auto& p : credible_regions_) {
        if (p.second.somatic) {
            using utils::to_string;
            if (!map_vafs_.empty()) {
                record.set_format(p.first, "MAP_VAF", to_string(map_vafs_.at(p.first)));
            }
            record.set_format(p.first, "VAF_CR", {to_string(p.second.somatic->first), to_string(p.second.somatic->second)});
        } else {
            if (!map_vafs_.empty()) {
                record.set_format_missing(p.first, "MAP_VAF");
            }
            record.set_format_missing(p.first, "VAF_CR");
        }
    }
}

std::unique_ptr<Call> SomaticCall::do_clone() const
{
    return std::make_unique<SomaticCall>(*this);
}

} // namespace octopus
