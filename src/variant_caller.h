//
//  variant_caller.h
//  Octopus
//
//  Created by Daniel Cooke on 29/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_caller__
#define __Octopus__variant_caller__

#include <vector>
#include <unordered_map>
#include <map>
#include <functional>

#include "allele.h"
#include "variant.h"
#include "variational_bayes_genotype_model.h"
#include "haplotype_phaser.h"

template <typename SampleIdType, typename RealType=double>
using VariantCalls = std::unordered_map<SampleIdType, std::map<Variant, RealType>>;

VariantCalls<HaplotypePhaser::SampleIdType, HaplotypePhaser::RealType>
call_variants(const std::vector<HaplotypePhaser::SampleIdType>& the_samples,
              const HaplotypePhaser::PhasedRegions& the_phased_regions,
              const std::vector<Variant>& the_candidates, const VariationalBayesGenotypeModel& the_model);

#endif /* defined(__Octopus__variant_caller__) */
