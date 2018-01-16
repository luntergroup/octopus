// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_realigner_hpp
#define read_realigner_hpp

#include <vector>

#include "basics/aligned_read.hpp"
#include "basics/cigar_string.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/haplotype_likelihood_model.hpp"

namespace octopus {

Haplotype expand_for_realignment(const Haplotype& haplotype, const std::vector<AlignedRead>& reads);

std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
                                 HaplotypeLikelihoodModel model);

std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype);

} // namespace

#endif
