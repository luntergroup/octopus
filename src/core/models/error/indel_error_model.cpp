// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "indel_error_model.hpp"

#include "core/types/haplotype.hpp"

namespace octopus {

std::unique_ptr<IndelErrorModel> IndelErrorModel::clone() const
{
    return do_clone();
}

IndelErrorModel::PenaltyType
IndelErrorModel::evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalities) const
{
    return do_evaluate(haplotype, gap_open_penalities);
}

} // namespace octopus
