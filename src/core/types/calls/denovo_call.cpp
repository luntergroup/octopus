// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_call.hpp"

namespace octopus {

void DenovoCall::decorate(VcfRecord::Builder& record) const
{
    record.set_denovo();
}

std::unique_ptr<Call> DenovoCall::do_clone() const
{
    return std::make_unique<DenovoCall>(*this);
}

} // namespace octopus
