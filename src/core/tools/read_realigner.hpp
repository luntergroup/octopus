// Copyright (c) 2015-2018 Daniel Cooke
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

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype);
void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model);
std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
                                 HaplotypeLikelihoodModel model);
std::vector<AlignedRead> realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype);

void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype);
std::vector<AlignedRead> safe_realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype);

CigarString rebase(const CigarString& read_to_haplotype, const CigarString& haplotype_to_reference);
void rebase(std::vector<AlignedRead>& reads, const Haplotype& haplotype);

void realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype);
std::vector<AlignedRead> realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype);

void safe_realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype);
std::vector<AlignedRead> safe_realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype);

} // namespace

#endif
