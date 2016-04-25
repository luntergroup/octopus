//
//  variant_call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_call.hpp"

#include <algorithm>
#include <unordered_map>
#include <utility>

#include "mappable.hpp"
#include "reference_genome.hpp"
#include "genotype.hpp"
#include "variant.hpp"

namespace Octopus
{
    const GenomicRegion& VariantCall::get_region() const noexcept
    {
        return mapped_region(variant_);
    }
    
    const Allele& VariantCall::get_reference() const noexcept
    {
        return variant_.get_ref_allele();
    }
    
    const Allele& VariantCall::get_alternative() const noexcept
    {
        return variant_.get_alt_allele();
    }
    
    void VariantCall::parsimonise(const ReferenceGenome& reference)
    {
        if (is_parsimonious(variant_)) return;
        
        if (all_genotypes_are_self_contained()) {
            auto parsimonised_variant = make_parsimonious(variant_, reference);
            
            const std::unordered_map<Allele, Allele> parsimonised_alleles {
                {variant_.get_ref_allele(), parsimonised_variant.get_ref_allele()},
                {variant_.get_alt_allele(), parsimonised_variant.get_alt_allele()}
            };
            
            variant_ = std::move(parsimonised_variant);
            
            for (auto& p : genotype_calls_) {
                Genotype<Allele>& genotype {p.second.genotype};
                
                Genotype<Allele> parsimonised_genotype {genotype.ploidy()};
                
                for (const Allele& allele : genotype) {
                    parsimonised_genotype.emplace(parsimonised_alleles.at(allele));
                }
                
                genotype = std::move(parsimonised_genotype);
            }
        } else {
            
        }
    }
    
    bool is_in(const Allele& allele, const Variant& variant)
    {
        return allele == variant.get_ref_allele() || allele == variant.get_alt_allele();
    }
    
    bool contains(const Genotype<Allele>& genotype, const Variant& variant)
    {
        return std::all_of(std::cbegin(genotype), std::cend(genotype),
                           [&variant] (const Allele& allele) {
                               return is_in(allele, variant);
                           });
    }
    
    bool VariantCall::all_genotypes_are_self_contained() const
    {
        return std::all_of(std::cbegin(genotype_calls_), std::cend(genotype_calls_),
                           [this] (const auto& p) {
                               return contains(p.second.genotype, variant_);
                           });
    }
} // namespace Octopus
