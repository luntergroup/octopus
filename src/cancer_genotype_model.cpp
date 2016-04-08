//
//  cancer_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_genotype_model.hpp"

#include <array>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <cstddef>
#include <memory>
#include <cmath>
#include <cassert>
#include <iostream>

#include <boost/math/special_functions/digamma.hpp>

#include "maths.hpp"
#include "logging.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
    // public methods
    
    Cancer::Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                   const unsigned ploidy, Priors priors)
    :
    Cancer {std::move(samples), normal_sample, ploidy, std::move(priors), AlgorithmParameters {}}
    {}
    
    Cancer::Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                   const unsigned ploidy, Priors priors, AlgorithmParameters parameters)
    :
    samples_ {std::move(samples)},
    normal_sample_ {normal_sample},
    ploidy_ {ploidy},
    priors_ {std::move(priors)},
    parameters_ {parameters}
    {
        assert(ploidy_ > 0);
        const auto it = std::find(std::begin(samples_), std::end(samples_), normal_sample_);
        assert(it != std::end(samples_));
        std::iter_swap(std::begin(samples_), it);
    }
    
    using ProbabilityVector    = std::vector<double>;
    using LogProbabilityVector = std::vector<double>;
    
    struct VariationalBayesParameters
    {
        explicit VariationalBayesParameters(double epsilon, unsigned max_iterations);
        double epsilon;
        unsigned max_iterations;
    };
    
    // Key non-member method declarations
    
    LogProbabilityVector calculate_log_priors(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                              const CoalescentModel& germline_prior_model);
    
    template <unsigned K>
    Cancer::InferredLatents
    run_variational_bayes(const std::vector<SampleIdType>& samples,
                          std::vector<CancerGenotype<Haplotype>>&& genotypes,
                          const Cancer::Priors& priors,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods,
                          const VariationalBayesParameters& params);
    
    // Cancer public
    
    Cancer::InferredLatents
    Cancer::infer_latents(const std::vector<Haplotype>& haplotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        assert(!haplotypes.empty());
        
        auto genotypes = generate_all_cancer_genotypes(haplotypes, ploidy_);
        
        assert(!genotypes.empty());
        
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "There are " << genotypes.size() << " initial cancer genotypes";
        }
        
        VariationalBayesParameters vb_params {parameters_.epsilon, parameters_.max_iterations};
        
        assert(ploidy_ < 3);
        
        if (ploidy_ == 1) {
            return run_variational_bayes<2>(samples_, std::move(genotypes), priors_,
                                            haplotype_likelihoods, vb_params);
        }
        
        return run_variational_bayes<3>(samples_, std::move(genotypes), priors_,
                                        haplotype_likelihoods, vb_params);
    }
    
    // Compressed types used by the Variational Bayes model
    
    template <unsigned K>
    using CompressedAlpha = std::array<double, K>;
    
    template <unsigned K>
    using CompressedAlphas = std::vector<CompressedAlpha<K>>;
    
    class ReadLikelihoods
    {
    public:
        using BaseType = HaplotypeLikelihoodCache::ReadProbabilities;
        
        ReadLikelihoods() = default;
        explicit ReadLikelihoods(const BaseType&);
        ~ReadLikelihoods() = default;
        
        ReadLikelihoods(const ReadLikelihoods&)            = default;
        ReadLikelihoods& operator=(const ReadLikelihoods&) = default;
        ReadLikelihoods(ReadLikelihoods&&)                 = default;
        ReadLikelihoods& operator=(ReadLikelihoods&&)      = default;
        
        void operator=(const BaseType&);
        void operator=(std::reference_wrapper<const BaseType>);
        
        std::size_t size() const noexcept;
        
        BaseType::const_iterator begin() const noexcept;
        BaseType::const_iterator end() const noexcept;
        
        double operator[](const std::size_t n) const noexcept;
        
    private:
        const BaseType* likelihoods;
    };
    
    template <unsigned K>
    using CompressedGenotype = std::array<ReadLikelihoods, K>;
    template <unsigned K>
    using CompressedGenotypes = std::vector<CompressedGenotype<K>>;
    template <unsigned K>
    using CompressedReadLikelihoods = std::vector<CompressedGenotypes<K>>;
    
    template <unsigned K>
    struct CompressedInferredLatents
    {
        ProbabilityVector genotype_posteriors;
        CompressedAlphas<K> alphas;
    };
    
    template <unsigned K>
    using Tau = std::array<double, K>;
    
    template <unsigned K>
    using ResponsabilityVector = std::vector<Tau<K>>;
    
    template <unsigned K>
    using ResponsabilityVectors = std::vector<ResponsabilityVector<K>>;
    
    // non-member methods
    
    void exp(ProbabilityVector& result, const LogProbabilityVector& log_probabilities)
    {
        std::transform(std::cbegin(log_probabilities), std::cend(log_probabilities),
                       std::begin(result), [] (const auto lp) { return std::exp(lp); });
    }
    
    ProbabilityVector exp(const LogProbabilityVector& log_probabilities)
    {
        ProbabilityVector result(log_probabilities.size());
        exp(result, log_probabilities);
        return result;
    }
    
    template <typename Container>
    void normalise_logs(Container& logs)
    {
        const auto norm = Maths::log_sum_exp(logs);
        for (auto& p : logs) p -= norm;
    }
    
    double calculate_log_prior(const CancerGenotype<Haplotype>& genotype,
                               const CoalescentModel& germline_prior_model)
    {
        const auto germline_log_prior = std::log(germline_prior_model.evaluate(genotype.get_germline_genotype()));
        double log_prob_cancer_hap_given_germline {0}; // TODO - uniform for now
        return germline_log_prior + log_prob_cancer_hap_given_germline;
    }
    
    LogProbabilityVector calculate_log_priors(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                              const CoalescentModel& germline_prior_model)
    {
        LogProbabilityVector result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&germline_prior_model] (const auto& genotype) {
                           return calculate_log_prior(genotype, germline_prior_model);
                       });
        
        normalise_logs(result);
        
        return result;
    }
    
    template <unsigned K>
    CompressedAlpha<K> compress(const Cancer::Priors::GenotypeMixturesDirichletAlphas& alpha)
    {
        CompressedAlpha<K> result;
        std::copy_n(std::cbegin(alpha), K, std::begin(result));
        return result;
    }
    
    template <unsigned K>
    CompressedAlphas<K> flatten_priors(const Cancer::Priors& priors, const std::vector<SampleIdType>& samples)
    {
        CompressedAlphas<K> result(samples.size());
        
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                       [&priors] (const auto& sample) {
                           return compress<K>(priors.alphas.at(sample));
                       });
        
        return result;
    }
    
    template <unsigned K>
    CompressedGenotype<K>
    compress(const CancerGenotype<Haplotype>& genotype, const SampleIdType& sample,
             const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        CompressedGenotype<K> result;
        
        const Genotype<Haplotype>& germline_genotype {genotype.get_germline_genotype()};
        
        assert(germline_genotype.ploidy() == (K - 1));
        
        std::transform(std::cbegin(germline_genotype), std::cend(germline_genotype),
                       std::begin(result),
                       [&sample, &haplotype_likelihoods] (const Haplotype& haplotype)
                            -> std::reference_wrapper<const ReadLikelihoods::BaseType> {
                           return std::cref(haplotype_likelihoods.log_likelihoods(sample, haplotype));
                       });
        
        result.back() = haplotype_likelihoods.log_likelihoods(sample, genotype.get_cancer_element());
        
        return result;
    }
    
    template <unsigned K>
    CompressedGenotypes<K>
    compress(const std::vector<CancerGenotype<Haplotype>>& genotypes, const SampleIdType& sample,
             const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        CompressedGenotypes<K> result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&sample, &haplotype_likelihoods] (const auto& genotype) {
                           return compress<K>(genotype, sample, haplotype_likelihoods);
                       });
        
        return result;
    }
    
    template <unsigned K>
    CompressedReadLikelihoods<K>
    compress(const std::vector<CancerGenotype<Haplotype>>& genotypes,
             const std::vector<SampleIdType>& samples,
             const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        CompressedReadLikelihoods<K> result {};
        result.reserve(samples.size());
        
        std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                       [&genotypes, &haplotype_likelihoods] (const auto& sample) {
                           return compress<K>(genotypes, sample, haplotype_likelihoods);
                       });
        
        return result;
    }
    
    auto sum(const CompressedAlpha<2>& alpha)
    {
        return alpha.front() + alpha.back();
    }
    
    auto sum(const CompressedAlpha<3>& alpha)
    {
        return alpha[0] + alpha[1] + alpha[2];
    }
    
    template <unsigned K>
    auto sum(const CompressedAlpha<K>& alpha)
    {
        return std::accumulate(std::cbegin(alpha), std::cend(alpha), 0.0);
    }
    
    auto digamma_diff(const double a, const double b)
    {
        using boost::math::digamma;
        return digamma(a) - digamma(b);
    }
    
    template <unsigned K>
    double expectation(const ProbabilityVector& distribution,
                       const CompressedGenotypes<K>& likelihoods,
                       const unsigned k, const std::size_t n)
    {
        return std::inner_product(std::cbegin(distribution), std::cend(distribution),
                                  std::cbegin(likelihoods), 0.0,
                                  std::plus<void> {},
                                  [k, n] (const auto p, const auto& haplotype_likelihoods) {
                                      return p * haplotype_likelihoods[k][n];
                                  });
    }
    
    template <unsigned K>
    ResponsabilityVectors<K>
    init_responsabilities(const CompressedAlphas<K>& prior_alphas,
                          const ProbabilityVector& genotype_priors,
                          const CompressedReadLikelihoods<K>& read_likelihoods)
    {
        assert(!read_likelihoods.empty());
        assert(prior_alphas.size() == read_likelihoods.size());
        
        // notation follows documentation
        const auto S = read_likelihoods.size(); // num samples
        
        ResponsabilityVectors<K> result {};
        result.reserve(S);
        
        for (std::size_t s {0}; s < S; ++s) {
            std::array<double, K> al; // no need to keep recomputing this
            
            const auto a0 = sum(prior_alphas[s]);
            
            for (unsigned k {0}; k < K; ++k) {
                al[k] = digamma_diff(prior_alphas[s][k], a0);
            }
            
            const auto N = read_likelihoods[s][0][0].size(); // num reads in sample s
            
            ResponsabilityVector<K> read_responsabilities(N);
            
            std::array<double, K> ln_rho;
            
            for (std::size_t n {0}; n < N; ++n) {
                for (unsigned k {0}; k < K; ++k) {
                    ln_rho[k] = al[k] + expectation<K>(genotype_priors, read_likelihoods[s], k, n);
                }
                
                const auto ln_rho_norm = Maths::log_sum_exp(ln_rho);
                
                for (unsigned k {0}; k < K; ++k) {
                    read_responsabilities[n][k] = std::exp(ln_rho[k] - ln_rho_norm);
                }
            }
            
            result.emplace_back(std::move(read_responsabilities));
        }
        
        return result;
    }
    
    // same as init_responsabilities but in-place
    template <unsigned K>
    void update_responsabilities(ResponsabilityVectors<K>& result,
                                 const CompressedAlphas<K>& posterior_alphas,
                                 const ProbabilityVector& genotype_posteriors,
                                 const CompressedReadLikelihoods<K>& read_likelihoods)
    {
        const auto S = read_likelihoods.size();
        for (std::size_t s {0}; s < S; ++s) {
            std::array<double, K> al;
            
            const auto a0 = sum(posterior_alphas[s]);
            
            for (unsigned k {0}; k < K; ++k) {
                al[k] = digamma_diff(posterior_alphas[s][k], a0);
            }
            
            const auto N = read_likelihoods[s][0][0].size();
            
            std::array<double, K> ln_rho;
            
            for (std::size_t n {0}; n < N; ++n) {
                for (unsigned k {0}; k < K; ++k) {
                    ln_rho[k] = al[k] + expectation<K>(genotype_posteriors, read_likelihoods[s], k, n);
                }
                
                const auto ln_rho_norm = Maths::log_sum_exp(ln_rho);
                
                for (unsigned k {0}; k < K; ++k) {
                    result[s][n][k] = std::exp(ln_rho[k] - ln_rho_norm);
                }
            }
        }
    }
    
    template <unsigned K>
    double sum(const ResponsabilityVector<K>& log_taus, const unsigned k)
    {
        return std::accumulate(std::cbegin(log_taus), std::cend(log_taus), 0.0,
                               [k] (const auto curr, const auto& tau) {
                                   return curr + tau[k];
                               });
    }
    
    template <unsigned K>
    void update_alpha(CompressedAlpha<K>& alpha,
                      const CompressedAlpha<K>& prior_alpha,
                      const ResponsabilityVector<K>& log_taus)
    {
        for (unsigned k {0}; k < K; ++k) {
            alpha[k] = prior_alpha[k] + sum<K>(log_taus, k);
        }
    }
    
    template <unsigned K>
    void update_alphas(CompressedAlphas<K>& alphas,
                       const CompressedAlphas<K>& prior_alphas,
                       const ResponsabilityVectors<K>& responsabilities)
    {
        const auto S = alphas.size();
        assert(S == prior_alphas.size() && S == responsabilities.size());
        for (std::size_t s {0}; s < S; ++s) {
            update_alpha<K>(alphas[s], prior_alphas[s], responsabilities[s]);
        }
    }
    
    template <unsigned K>
    double marginalise(const ResponsabilityVectors<K>& responsabilities,
                       const CompressedReadLikelihoods<K>& read_likelihoods,
                       const std::size_t g)
    {
        double result {0};
        
        const auto S = read_likelihoods.size(); // num samples
        
        assert(S == responsabilities.size());
        
        for (std::size_t s {0}; s < S; ++s) {
            const auto N = read_likelihoods[s][0][0].size(); // num reads in sample s
            
            for (unsigned k {0}; k < K; ++k) {
                for (std::size_t n {0}; n < N; ++n) {
                    result += responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
                }
            }
        }
        
        return result;
    }
    
    template <unsigned K>
    void update_genotype_log_posteriors(LogProbabilityVector& result,
                                        const LogProbabilityVector& genotype_log_priors,
                                        const ResponsabilityVectors<K>& responsabilities,
                                        const CompressedReadLikelihoods<K>& read_likelihoods)
    {
        for (std::size_t g {0}; g < result.size(); ++g) {
            result[g] = genotype_log_priors[g] + marginalise<K>(responsabilities, read_likelihoods, g);
        }
        normalise_logs(result);
    }
    
    auto max_change(const CompressedAlpha<2>& lhs, const CompressedAlpha<2>& rhs)
    {
        return std::max(std::abs(lhs.front() - lhs.front()), std::abs(lhs.back() - lhs.back()));
    }
    
    auto max_change(const CompressedAlpha<3>& lhs, const CompressedAlpha<3>& rhs)
    {
        return std::max({std::abs(lhs[0] - lhs[0]), std::abs(lhs[1] - lhs[1]), std::abs(lhs[2] - lhs[2])});
    }
    
    template <unsigned K>
    auto max_change(const CompressedAlpha<K>& lhs, const CompressedAlpha<K>& rhs)
    {
        double result {0};
        
        for (std::size_t k {0}; k < K; ++k) {
            const auto curr = std::abs(lhs[k] - rhs[k]);
            if (curr > result) result = curr;
        }
        
        return result;
    }
    
    template <unsigned K>
    auto max_change(const CompressedAlphas<K>& prior_alphas,
                    const CompressedAlphas<K>& posterior_alphas)
    {
        assert(prior_alphas.size() == posterior_alphas.size());
        
        const auto S = prior_alphas.size();
        
        double result {0};
        
        for (std::size_t s {0}; s < S; ++s) {
            const auto curr = max_change<K>(prior_alphas[s], posterior_alphas[s]);
            if (curr > result) result = curr;
        }
        
        return result;
    }
    
    template <unsigned K>
    std::pair<bool, double> check_convergence(const CompressedAlphas<K>& prior_alphas,
                                              const CompressedAlphas<K>& posterior_alphas,
                                              const double prev_max_change,
                                              const double epsilon)
    {
        const auto new_max_change = max_change<K>(prior_alphas, posterior_alphas);
        return std::make_pair(std::abs(new_max_change - prev_max_change) < epsilon, new_max_change);
    }
    
    // The main algorithm
    template <unsigned K>
    CompressedInferredLatents<K>
    run_variational_bayes(const CompressedAlphas<K>& prior_alphas,
                          const LogProbabilityVector& genotype_log_priors,
                          const CompressedReadLikelihoods<K>& read_likelihoods,
                          const VariationalBayesParameters params)
    {
        assert(!prior_alphas.empty());
        assert(!genotype_log_priors.empty());
        assert(!read_likelihoods.empty());
        assert(prior_alphas.size() == read_likelihoods.size()); // num samples
        assert(read_likelihoods.front().size() == genotype_log_priors.size()); // num genotypes
        assert(params.max_iterations > 0);
        
        // initialise everything first
        auto genotype_log_posteriors = genotype_log_priors;
        auto genotype_posteriors     = exp(genotype_log_posteriors);
        
        auto posterior_alphas = prior_alphas;
        
        auto responsabilities = init_responsabilities<K>(posterior_alphas, genotype_posteriors,
                                                         read_likelihoods);
        
        assert(responsabilities.size() == read_likelihoods.size()); // num samples
        assert(!responsabilities.front().empty());
        
        
        bool is_converged;
        double max_change;
        // main loop
        for (unsigned i {0}; i < params.max_iterations; ++i) {
            std::cout << "VB Iteration " << i << std::endl;
            update_genotype_log_posteriors<K>(genotype_log_posteriors, genotype_log_priors,
                                              responsabilities, read_likelihoods);
            
            exp(genotype_posteriors, genotype_log_posteriors);
            
            update_alphas<K>(posterior_alphas, prior_alphas, responsabilities);
            
            update_responsabilities<K>(responsabilities, posterior_alphas, genotype_posteriors,
                                       read_likelihoods);
            
            std::tie(is_converged, max_change) = check_convergence<K>(prior_alphas, posterior_alphas,
                                                                      max_change, params.epsilon);
            
            if (is_converged) {
                break;
            }
        }
        
        return {std::move(genotype_posteriors), std::move(posterior_alphas)};
    }
    
    Cancer::InferredLatents::GenotypePosteriorMap
    expand(std::vector<CancerGenotype<Haplotype>>&& genotypes,
           LogProbabilityVector&& genotype_log_posteriors)
    {
        Cancer::InferredLatents::GenotypePosteriorMap result {};
        
        std::transform(std::make_move_iterator(std::begin(genotypes)),
                       std::make_move_iterator(std::end(genotypes)),
                       std::begin(genotype_log_posteriors),
                       std::inserter(result, std::begin(result)),
                       [] (auto&& g, auto p) {
                           return std::make_pair(std::move(g), p);
                       });
        
        return result;
    }
    
    template <unsigned K>
    Cancer::InferredLatents::GenotypeMixturesDirichletAlphas expand(CompressedAlpha<K>& alpha)
    {
        return Cancer::InferredLatents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
    }
    
    template <unsigned K>
    Cancer::InferredLatents::GenotypeMixturesDirichletAlphaMap
    expand(const std::vector<SampleIdType>& samples, CompressedAlphas<K>&& alphas)
    {
        Cancer::InferredLatents::GenotypeMixturesDirichletAlphaMap result {};
        
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                       std::inserter(result, std::begin(result)),
                       [] (const auto& sample, auto&& compressed_alpha) {
                           return std::make_pair(sample, expand<K>(compressed_alpha));
                       });
        
        return result;
    }
    
    template <unsigned K>
    Cancer::InferredLatents
    expand(const std::vector<SampleIdType>& samples,
           std::vector<CancerGenotype<Haplotype>>&& genotypes,
           CompressedInferredLatents<K>&& inferred_latents)
    {
        return Cancer::InferredLatents {
            expand(std::move(genotypes), std::move(inferred_latents.genotype_posteriors)),
            expand<K>(samples, std::move(inferred_latents.alphas))
        };
    }
    
    template <unsigned K>
    Cancer::InferredLatents
    run_variational_bayes(const std::vector<SampleIdType>& samples,
                          std::vector<CancerGenotype<Haplotype>>&& genotypes,
                          const Cancer::Priors& priors,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods,
                          const VariationalBayesParameters& params)
    {
        const auto prior_alphas = flatten_priors<K>(priors, samples);
        
        const auto genotype_log_priors = calculate_log_priors(genotypes, priors.genotype_prior_model);
        
        for (std::size_t g {0}; g < genotypes.size(); ++g) {
            ::debug::print_variant_alleles(genotypes[g]); std::cout << " " << genotype_log_priors[g] << '\n';
        }
        //exit(0);
        
        const auto read_likelihoods = compress<K>(genotypes, samples, haplotype_likelihoods);
        
        auto inferred_latents = run_variational_bayes<K>(prior_alphas, genotype_log_priors,
                                                         read_likelihoods, params);
        
        return expand<K>(samples, std::move(genotypes), std::move(inferred_latents));
    }
    
    VariationalBayesParameters::VariationalBayesParameters(const double epsilon,
                                                           const unsigned max_iterations)
    :
    epsilon {epsilon},
    max_iterations {max_iterations}
    {}
    
    ReadLikelihoods::ReadLikelihoods(const BaseType& underlying_likelihoods)
    :
    likelihoods {std::addressof(underlying_likelihoods)}
    {}
    
    void ReadLikelihoods::operator=(const BaseType& other)
    {
        likelihoods = std::addressof(other);
    }
    
    void ReadLikelihoods::operator=(std::reference_wrapper<const BaseType> other)
    {
        likelihoods = std::addressof(other.get());
    }
    
    std::size_t ReadLikelihoods::size() const noexcept
    {
        return likelihoods->size();
    }
    
    ReadLikelihoods::BaseType::const_iterator ReadLikelihoods::begin() const noexcept
    {
        return likelihoods->begin();
    }
    
    ReadLikelihoods::BaseType::const_iterator ReadLikelihoods::end() const noexcept
    {
        return likelihoods->end();
    }
    
    double ReadLikelihoods::operator[](const std::size_t n) const noexcept
    {
        return likelihoods->operator[](n);
    }
    
    namespace debug
    {
        
    } // namespace debug
    
    } // namespace GenotypeModel
} // namespace Octopus
