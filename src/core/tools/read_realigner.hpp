// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_realigner_hpp
#define read_realigner_hpp

#include <vector>

#include <boost/optional.hpp>

#include "basics/aligned_read.hpp"
#include "basics/cigar_string.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "utils/thread_pool.hpp"

namespace octopus {

Haplotype expand_for_realignment(const Haplotype& haplotype, const std::vector<AlignedRead>& reads,
                                 const HaplotypeLikelihoodModel& model);
Haplotype expand_for_realignment(const Haplotype& haplotype, const std::vector<AlignedRead>& reads);

void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype,
             HaplotypeLikelihoodModel model, boost::optional<ThreadPool&> workers = boost::none);
void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype,
             boost::optional<ThreadPool&> workers = boost::none);
void realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype,
             HaplotypeLikelihoodModel model,
             std::vector<HaplotypeLikelihoodModel::LogProbability>& log_likelihoods,
             boost::optional<ThreadPool&> workers = boost::none);

std::vector<AlignedRead>
realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model,
        boost::optional<ThreadPool&> workers = boost::none);
std::vector<AlignedRead>
realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers = boost::none);

std::vector<AlignedRead>
realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype,
        HaplotypeLikelihoodModel model,
        std::vector<HaplotypeLikelihoodModel::LogProbability>& log_likelihoods,
        boost::optional<ThreadPool&> workers = boost::none);

void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype,
                  HaplotypeLikelihoodModel model, 
                  boost::optional<ThreadPool&> workers = boost::none);
void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers = boost::none);
void safe_realign(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model,
                  std::vector<HaplotypeLikelihoodModel::LogProbability>& log_likelihoods,
                  boost::optional<ThreadPool&> workers = boost::none);

std::vector<AlignedRead>
safe_realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model,
             boost::optional<ThreadPool&> workers = boost::none);
std::vector<AlignedRead>
safe_realign(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers = boost::none);

CigarString rebase(const CigarString& read_to_haplotype, const CigarString& haplotype_to_reference);
void rebase(std::vector<AlignedRead>& reads, const Haplotype& haplotype, boost::optional<ThreadPool&> workers = boost::none);

void realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model);
void realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype);

std::vector<AlignedRead>
realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model);
std::vector<AlignedRead>
realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype);

void safe_realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model);
void safe_realign_to_reference(std::vector<AlignedRead>& reads, const Haplotype& haplotype);

std::vector<AlignedRead>
safe_realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype, HaplotypeLikelihoodModel model);
std::vector<AlignedRead>
safe_realign_to_reference(const std::vector<AlignedRead>& reads, const Haplotype& haplotype);

} // namespace

#endif
