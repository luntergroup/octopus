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
    
    struct DummyGenerator
    {
        DummyGenerator() = delete;
        DummyGenerator(const char dummy) : dummy_ {dummy} {}
        
        char operator()(const GenomicRegion& region) const
        {
            return dummy_;
        }
    private:
        char dummy_;
    };
    
    void VariantCall::parsimonise(const char dummy_base)
    {
        if (is_parsimonious(variant_)) return;
        
        auto parsimonised_variant = make_parsimonious(variant_, DummyGenerator {dummy_base});
        
        const std::unordered_map<Allele, Allele> parsimonised_alleles {
            {variant_.get_ref_allele(), parsimonised_variant.get_ref_allele()},
            {variant_.get_alt_allele(), parsimonised_variant.get_alt_allele()}
        };
        
        const auto has_variant_shifted = begins_before(parsimonised_variant, variant_);
        
        variant_ = std::move(parsimonised_variant);
        
        for (auto& p : genotype_calls_) {
            Genotype<Allele>& genotype {p.second.genotype};
            
            Genotype<Allele> parsimonised_genotype {genotype.ploidy()};
            
            for (const Allele& allele : genotype) {
                if (parsimonised_alleles.count(allele) == 1) {
                    parsimonised_genotype.emplace(parsimonised_alleles.at(allele));
                } else {
                    if (has_variant_shifted) {
                        auto old_sequence = allele.get_sequence();
                        old_sequence.insert(std::begin(old_sequence), dummy_base);
                        Allele new_allele {mapped_region(variant_), std::move(old_sequence)};
                        parsimonised_genotype.emplace(std::move(new_allele));
                    } else {
                        parsimonised_genotype.emplace(allele);
                    }
                }
            }
            
            genotype = std::move(parsimonised_genotype);
        }
    }
    
    void VariantCall::parsimonise(const ReferenceGenome& reference)
    {
        if (is_parsimonious(variant_)) return;
        
        auto parsimonised_variant = make_parsimonious(variant_, reference);
        
        const std::unordered_map<Allele, Allele> parsimonised_alleles {
            {variant_.get_ref_allele(), parsimonised_variant.get_ref_allele()},
            {variant_.get_alt_allele(), parsimonised_variant.get_alt_allele()}
        };
        
        const auto has_variant_shifted = begins_before(parsimonised_variant, variant_);
        
        char reference_base;
        if (has_variant_shifted) {
            reference_base = reference.get_sequence(head_position(parsimonised_variant)).front();
        }
        
        variant_ = std::move(parsimonised_variant);
        
        for (auto& p : genotype_calls_) {
            Genotype<Allele>& genotype {p.second.genotype};
            
            Genotype<Allele> parsimonised_genotype {genotype.ploidy()};
            
            for (const Allele& allele : genotype) {
                if (parsimonised_alleles.count(allele) == 1) {
                    parsimonised_genotype.emplace(parsimonised_alleles.at(allele));
                } else {
                    if (has_variant_shifted) {
                        auto old_sequence = allele.get_sequence();
                        old_sequence.insert(std::begin(old_sequence), reference_base);
                        Allele new_allele {mapped_region(variant_), std::move(old_sequence)};
                        parsimonised_genotype.emplace(std::move(new_allele));
                    } else {
                        parsimonised_genotype.emplace(allele);
                    }
                }
            }
            
            genotype = std::move(parsimonised_genotype);
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
