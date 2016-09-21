// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_call.hpp"

namespace octopus {

void DenovoCall::decorate(VcfRecord::Builder& record) const
{
    record.set_denovo();
}

} // namespace octopus
