// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_call.hpp"

#include "utils/string_utils.hpp"

namespace octopus {

void DenovoCall::decorate(VcfRecord::Builder& record) const
{
    record.set_denovo();
    if (posterior_) {
        record.set_info("PP", utils::to_string(posterior_->score()));
    }
}

std::unique_ptr<Call> DenovoCall::do_clone() const
{
    return std::make_unique<DenovoCall>(*this);
}

} // namespace octopus
