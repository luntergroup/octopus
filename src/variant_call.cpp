//
//  variant_call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_call.hpp"

#include <algorithm>
#include <iterator>
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
    
    void VariantCall::replace_called_alleles(const char old_base, const char replacement_base)
    {
        const auto& ref_seq = ref_sequence(variant_);
        const auto& alt_seq = alt_sequence(variant_);
        
        const auto it1 = std::find(std::cbegin(ref_seq), std::cend(ref_seq), old_base);
        const auto it2 = std::find(std::cbegin(alt_seq), std::cend(alt_seq), old_base);
        
        if (it1 == std::cend(ref_seq) && it2 == std::cend(alt_seq)) {
            return;
        }
        
        Allele new_ref, new_alt;
        
        if (it1 != std::cend(ref_seq)) {
            Allele::SequenceType new_sequence {};
            new_sequence.reserve(ref_seq.size());
            
            new_sequence.insert(std::end(new_sequence), std::cbegin(ref_seq), it1);
            
            std::replace_copy(it1, std::cend(ref_seq), std::back_inserter(new_sequence),
                              old_base, replacement_base);
            
            new_ref = Allele {variant_.get_region(), std::move(new_sequence)};
        } else {
            new_ref = variant_.get_ref_allele();
        }
        
        if (it2 != std::cend(alt_seq)) {
            Allele::SequenceType new_sequence {};
            new_sequence.reserve(alt_seq.size());
            
            new_sequence.insert(std::end(new_sequence), std::cbegin(alt_seq), it2);
            
            std::replace_copy(it2, std::cend(alt_seq), std::back_inserter(new_sequence),
                              old_base, replacement_base);
            
            new_alt = Allele {variant_.get_region(), std::move(new_sequence)};
        } else {
            new_alt = variant_.get_alt_allele();
        }
        
        variant_ = Variant {std::move(new_ref), std::move(new_alt)};
    }
    
    void VariantCall::replace(const Allele& old, Allele replacement)
    {
        if (variant_.get_ref_allele() == old) {
            variant_ = Variant {std::move(replacement), variant_.get_alt_allele()};
        } else if (variant_.get_alt_allele() == old) {
            variant_ = Variant {variant_.get_ref_allele(), std::move(replacement)};
        }
    }
    
    bool matches(const Allele::SequenceType& lhs, const Allele::SequenceType& rhs,
                 const char ignoring)
    {
        auto it = std::mismatch(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs),
                                [ignoring] (const auto a, const auto b) {
                                    return a == b || a == ignoring || b == ignoring;
                                });
        
        return it.first == std::cend(lhs);
    }
    
    bool matches(const Allele& allele, const Variant& variant, const char ignoring)
    {
        if (allele == variant.get_ref_allele()
            || allele == variant.get_alt_allele()) return true;
        
        if (allele.get_sequence().size() == ref_sequence_size(variant)) {
            return matches(ref_sequence(variant), allele.get_sequence(), ignoring);
        }
        
        if (allele.get_sequence().size() == alt_sequence_size(variant)) {
            return matches(alt_sequence(variant), allele.get_sequence(), ignoring);
        }
        
        return false;
    }
    
    void VariantCall::replace_uncalled_genotype_alleles(const Allele& replacement,
                                                        const char ignoring)
    {
        for (auto& p : genotype_calls_) {
            auto it = std::find_if_not(std::cbegin(p.second.genotype), std::cend(p.second.genotype),
                                       [this, ignoring] (const Allele& allele) {
                                           return matches(allele, variant_, ignoring);
                                       });
            
            if (it != std::cend(p.second.genotype)) {
                Genotype<Allele> new_genotype {p.second.genotype.ploidy()};
                
                std::for_each(std::cbegin(p.second.genotype), it,
                              [&new_genotype] (const Allele& allele) {
                                  new_genotype.emplace(allele);
                              });
                new_genotype.emplace(replacement);
                std::for_each(std::next(it), std::cend(p.second.genotype),
                              [this, &new_genotype, &replacement, ignoring] (const Allele& allele) {
                                  if (matches(allele, variant_, ignoring)) {
                                      new_genotype.emplace(allele);
                                  } else {
                                      new_genotype.emplace(replacement);
                                  }
                              });
                
                p.second.genotype = std::move(new_genotype);
            }
        }
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
    
    bool VariantCall::parsimonise(const char dummy_base)
    {
        if (is_parsimonious(variant_)) return false;
        
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
        
        return has_variant_shifted;
    }
    
    bool VariantCall::parsimonise(const ReferenceGenome& reference)
    {
        if (is_parsimonious(variant_)) return false;
        
        auto parsimonised_variant = make_parsimonious(variant_, reference);
        
        const std::unordered_map<Allele, Allele> parsimonised_alleles {
            {variant_.get_ref_allele(), parsimonised_variant.get_ref_allele()},
            {variant_.get_alt_allele(), parsimonised_variant.get_alt_allele()}
        };
        
        const auto has_variant_shifted = begins_before(parsimonised_variant, variant_);
        
        char reference_base {};
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
        
        return has_variant_shifted;
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
