// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "germline_variant_call.hpp"

namespace octopus {

void GermlineVariantCall::decorate(VcfRecord::Builder& record) const
{
    
}

std::unique_ptr<Call> GermlineVariantCall::do_clone() const
{
    return std::make_unique<GermlineVariantCall>(*this);
}

} // namespace octopus
