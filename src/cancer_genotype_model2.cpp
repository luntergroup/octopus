//
//  cancer_genotype_model2.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

//#include "cancer_genotype_model.hpp"
//
//#include <array>
//#include <numeric>
//#include <algorithm>
//#include <iterator>
//#include <cmath>
//
//#include "common.hpp"
//#include "maths.hpp"
//#include "read_utils.hpp"
//
//#include <iostream>   // DEBUG
//
//namespace Octopus
//{
//    namespace GenotypeModel
//    {
//        // public methods
//        
//        Cancer::Cancer(SampleIdType normal_sample, unsigned max_em_iterations, double em_epsilon)
//        :
//        max_em_iterations_ {max_em_iterations},
//        em_epsilon_ {em_epsilon},
//        normal_sample_ {std::move(normal_sample)}
//        {}
//        
//        // non members
//        
//        using HaplotypeLogFrequencies           = std::unordered_map<Haplotype, double>;
//        
//        struct SampleGenotypeLogWeightsCounts
//        {
//            SampleGenotypeLogWeightsCounts() = default;
//            
//            std::array<double, 3> log_counts;
//            double log_count_sum;
//        };
//        
//        using SampleGenotypeLogWeightsCounts    = std::array<double, 3>;
//        using GenotypeWeightsLogCounts          = std::unordered_map<SampleIdType, SampleGenotypeWeightsCounts>;
//        
//        using SampleGenotypeLogWeights          = std::array<double, 3>;
//        using GenotypeLogWeights                = std::unordered_map<SampleIdType, SampleGenotypeWeights>;
//        
//        using GenotypeWeightLogResponsibilities = std::unordered_map<SampleIdType, std::vector<std::array<double, 3>>>;
//        
//        struct GenotypeLogPosterior
//        {
//            GenotypeLogPosterior() = delete;
//            GenotypeLogPosterior(const CancerGenotype<Haplotype>& genotype, double log_posterior)
//            : genotype {genotype}, log_posterior {log_posterior} {}
//            
//            const CancerGenotype<Haplotype>& genotype;
//            double log_posterior;
//        };
//        
//        using GenotypeLogPosteriors = std::vector<GenotypeLogPosterior>;
//        
//        GenotypeLogWeights init_genotype_log_weights(const GenotypeWeightsLogCounts& genotype_weight_log_counts)
//        {
//            GenotypeLogWeights result {};
//            result.reserve(weight_counts.size());
//            
//            for (const auto& sample_weight_log_counts : genotype_weight_log_counts) {
//                const auto& sample_priors = sample_weight_log_counts.second;
//                
//                auto sample_result = sample_priors.log_counts;
//                
//                sample_result[0] -= sample_priors.log_count_sum;
//                sample_result[1] -= sample_priors.log_count_sum;
//                sample_result[2] -= sample_priors.log_count_sum;
//                
//                result.emplace(sample_weight_counts.first, sample_result);
//            }
//            
//            return result;
//        }
//        
//        double
//        update_genotype_log_weights(GenotypeLogWeights& current_genotype_log_weights,
//                                    const GenotypeWeightsLogCounts& genotype_weight_log_priors,
//                                    const GenotypeWeightLogResponsibilities& genotype_weight_log_responsibilities)
//        {
//            double max_log_weight_change {0};
//            
//            for (const auto& current_log_weights : current_genotype_log_weights) {
//                const auto& sample = current_log_weights.first;
//                const auto& z = genotype_weight_responsibilities.at(sample);
//                
//                const auto& sample_weight_log_priors = genotype_weight_log_counts.at(sample);
//                
//                auto norm = Maths::log_sum_exp(z.size(), sample_weight_log_priors.log_sum);
//                
//                for (unsigned k {}; k < 3; ++k) {
//                    curr[k] = std::accumulate(std::cbegin(z), std::cend(z), 0.0,
//                                              [k] (double v, const auto& a) { return v + a[k]; });
//                    curr[k] += sample_prior_counts.second[k];
//                    curr[k] /= norm;
//                }
//                
//                result.emplace(sample, curr);
//            }
//            
//            return max_log_weight_change;
//        }
//        
//    } // namespace GenotypeModel
//} // namespace Octopus
