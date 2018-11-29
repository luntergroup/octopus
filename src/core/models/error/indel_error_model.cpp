// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "indel_error_model.hpp"

#include "core/types/haplotype.hpp"

namespace octopus {

std::unique_ptr<IndelErrorModel> IndelErrorModel::clone() const
{
    return do_clone();
}

void IndelErrorModel::set_penalties(const Haplotype& haplotype, PenaltyVector& gap_open_penalities, PenaltyType& gap_extend_penalty) const
{
    do_set_penalties(haplotype, gap_open_penalities, gap_extend_penalty);
}

} // namespace octopus
