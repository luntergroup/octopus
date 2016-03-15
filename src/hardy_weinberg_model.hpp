//
//  hardy_weinberg_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef hardy_weinberg_model_hpp
#define hardy_weinberg_model_hpp

#include <vector>
#include <cmath>

#include "maths.hpp"

namespace Octopus
{
    namespace detail
    {
        template <typename Genotype, typename Map>
        double log_hardy_weinberg_haploid(const Genotype& genotype,
                                          const Map& haplotype_frequencies)
        {
            return std::log(haplotype_frequencies.at(genotype[0]));
        }
        
        template <typename Genotype, typename Map>
        double log_hardy_weinberg_diploid(const Genotype& genotype,
                                          const Map& haplotype_frequencies)
        {
            if (genotype.is_homozygous()) {
                return 2 * std::log(haplotype_frequencies.at(genotype[0]));
            }
            
            static const double ln_2 {std::log(2.0)};
            
            return std::log(haplotype_frequencies.at(genotype[0]))
                        + std::log(haplotype_frequencies.at(genotype[1])) + ln_2;
        }
        
        template <typename Genotype, typename Map>
        double log_hardy_weinberg_triploid(const Genotype& genotype,
                                           const Map& haplotype_frequencies)
        {
            // TODO: optimise this case
            auto unique_haplotypes = genotype.copy_unique();
            
            std::vector<unsigned> occurences {};
            occurences.reserve(unique_haplotypes.size());
            
            double r {0};
            
            for (const auto& haplotype : unique_haplotypes) {
                auto num_occurences = genotype.count(haplotype);
                occurences.push_back(num_occurences);
                r += num_occurences * std::log(haplotype_frequencies.at(haplotype));
            }
            
            return Maths::log_multinomial_coefficient<double>(occurences) + r;
        }
        
        template <typename Genotype, typename Map>
        double log_hardy_weinberg_polyploid(const Genotype& genotype,
                                            const Map& haplotype_frequencies)
        {
            auto unique_haplotypes = genotype.copy_unique();
            
            std::vector<unsigned> occurences {};
            occurences.reserve(unique_haplotypes.size());
            
            double r {0};
            
            for (const auto& haplotype : unique_haplotypes) {
                auto num_occurences = genotype.count(haplotype);
                occurences.push_back(num_occurences);
                r += num_occurences * std::log(haplotype_frequencies.at(haplotype));
            }
            
            return Maths::log_multinomial_coefficient<double>(occurences) + r;
        }
    }
    
    // TODO: improve this, possible bottleneck in EM update at the moment
    template <typename Genotype, typename Map>
    double log_hardy_weinberg(const Genotype& genotype, const Map& haplotype_frequencies)
    {
        switch (genotype.ploidy()) {
            case 1 : return detail::log_hardy_weinberg_haploid(genotype, haplotype_frequencies);
            case 2 : return detail::log_hardy_weinberg_diploid(genotype, haplotype_frequencies);
            case 3 : return detail::log_hardy_weinberg_triploid(genotype, haplotype_frequencies);
            default: return detail::log_hardy_weinberg_polyploid(genotype, haplotype_frequencies);
        }
    }
}

#endif /* hardy_weinberg_model_hpp */
