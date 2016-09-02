// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
 
#include "trio_model.hpp"
 
#include <cassert>
 
#include "utils/maths.hpp"
#include "germline_likelihood_model.hpp"
 
namespace octopus { namespace model {
 
TrioModel::TrioModel(const Trio& trio,
              const CoalescentModel& genotype_prior_model,
              const DeNovoModel& mutation_model,
              boost::optional<logging::DebugLogger> debug_log)
: trio_ {trio}
, genotype_prior_model_ {genotype_prior_model}
, mutation_model_ {mutation_model}
, debug_log_ {debug_log}
{}
 
TrioModel::InferredLatents
TrioModel::evaluate(const GenotypeVector& maternal_genotypes, const GenotypeVector& paternal_genotypes,
                    const GenotypeVector& child_genotypes,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!maternal_genotypes.empty());
    assert(!paternal_genotypes.empty());
    assert(!child_genotypes.empty());
     
    const GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
     
     
     
    InferredLatents result {};
     
    return result;
}
 
} // namespace model
} // namespace octopus
