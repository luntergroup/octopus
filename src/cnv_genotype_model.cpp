//
//  cnv_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 12/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "cnv_genotype_model.hpp"

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
#include "individual_genotype_model.hpp"

namespace Octopus
{
namespace GenotypeModel
{
// public methods

CNV::CNV(std::vector<SampleIdType> samples, const unsigned ploidy, Priors priors)
:
CNV {std::move(samples), ploidy, std::move(priors), AlgorithmParameters {}}
{}

CNV::CNV(std::vector<SampleIdType> samples, const unsigned ploidy, Priors priors,
         AlgorithmParameters parameters)
:
samples_ {std::move(samples)},
ploidy_ {ploidy},
priors_ {std::move(priors)},
parameters_ {parameters}
{}

CNV::InferredLatents::InferredLatents(Latents&& posteriors, double approx_log_evidence)
:
posteriors {std::move(posteriors)},
approx_log_evidence {approx_log_evidence}
{}

namespace
{
    // Key non-member declarations
    
    using ProbabilityVector    = std::vector<double>;
    using LogProbabilityVector = std::vector<double>;
    
    struct VariationalBayesParameters
    {
        explicit VariationalBayesParameters(double epsilon, unsigned max_iterations);
        double epsilon;
        unsigned max_iterations;
    };
    
    template <std::size_t K>
    CNV::InferredLatents
    run_variational_bayes(const std::vector<SampleIdType>& samples,
                          std::vector<Genotype<Haplotype>>&& genotypes,
                          const CNV::Priors& priors,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods,
                          const VariationalBayesParameters& params);
}

// CNV public

CNV::InferredLatents
CNV::infer_latents(std::vector<Genotype<Haplotype>> genotypes,
                   const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    
    const VariationalBayesParameters vb_params {parameters_.epsilon, parameters_.max_iterations};
    
    assert(ploidy_ < 3);
    
    if (ploidy_ == 1) {
        return run_variational_bayes<1>(samples_, std::move(genotypes), priors_,
                                        haplotype_likelihoods, vb_params);
    }
    
    return run_variational_bayes<2>(samples_, std::move(genotypes), priors_,
                                    haplotype_likelihoods, vb_params);
}

namespace
{
    // Compressed types used by the Variational Bayes model
    
    template <std::size_t K>
    using CompressedAlpha = std::array<double, K>;
    
    template <std::size_t K>
    using CompressedAlphas = std::vector<CompressedAlpha<K>>;
    
    class ReadLikelihoods
    {
    public:
        using BaseType = HaplotypeLikelihoodCache::Likelihoods;
        
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
    
    template <std::size_t K>
    using CompressedGenotype = std::array<ReadLikelihoods, K>;
    template <std::size_t K>
    using CompressedGenotypes = std::vector<CompressedGenotype<K>>;
    template <std::size_t K>
    using CompressedReadLikelihoods = std::vector<CompressedGenotypes<K>>;
    
    template <std::size_t K>
    using Tau = std::array<double, K>;
    
    template <std::size_t K>
    using ResponsabilityVector = std::vector<Tau<K>>;
    
    template <std::size_t K>
    using ResponsabilityVectors = std::vector<ResponsabilityVector<K>>;
    
    template <std::size_t K>
    struct CompressedLatents
    {
        ProbabilityVector genotype_posteriors;
        LogProbabilityVector genotype_log_posteriors;
        CompressedAlphas<K> alphas;
        ResponsabilityVectors<K> responsabilities;
    };
    
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
    
    template <std::size_t K>
    CompressedAlpha<K> compress(const CNV::Priors::GenotypeMixturesDirichletAlphas& alpha)
    {
        CompressedAlpha<K> result;
        std::copy_n(std::cbegin(alpha), K, std::begin(result));
        return result;
    }
    
    template <std::size_t K>
    CompressedAlphas<K> flatten_priors(const CNV::Priors& priors,
                                       const std::vector<SampleIdType>& samples)
    {
        CompressedAlphas<K> result(samples.size());
        
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                       [&priors] (const auto& sample) {
                           return compress<K>(priors.alphas.at(sample));
                       });
        
        return result;
    }
    
    template <std::size_t K>
    CompressedGenotype<K>
    compress(const Genotype<Haplotype>& genotype, const SampleIdType& sample,
             const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        CompressedGenotype<K> result;
        
        std::transform(std::cbegin(genotype), std::cend(genotype),
                       std::begin(result),
                       [&sample, &haplotype_likelihoods] (const Haplotype& haplotype)
                       -> std::reference_wrapper<const ReadLikelihoods::BaseType> {
                           return std::cref(haplotype_likelihoods.log_likelihoods(sample, haplotype));
                       });
        
        return result;
    }
    
    template <std::size_t K>
    CompressedGenotypes<K>
    compress(const std::vector<Genotype<Haplotype>>& genotypes, const SampleIdType& sample,
             const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        CompressedGenotypes<K> result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&sample, &haplotype_likelihoods] (const auto& genotype) {
                           return compress<K>(genotype, sample, haplotype_likelihoods);
                       });
        
        return result;
    }
    
    template <std::size_t K>
    CompressedReadLikelihoods<K>
    compress(const std::vector<Genotype<Haplotype>>& genotypes,
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
    
    template <std::size_t K>
    auto sum(const CompressedAlpha<K>& alpha)
    {
        return std::accumulate(std::cbegin(alpha), std::cend(alpha), 0.0);
    }
    
    auto digamma_diff(const double a, const double b)
    {
        using boost::math::digamma;
        return digamma(a) - digamma(b);
    }
    
    template <std::size_t K>
    double expectation(const ProbabilityVector& distribution,
                       const CompressedGenotypes<K>& likelihoods,
                       const unsigned k, const std::size_t n)
    {
        return std::inner_product(std::cbegin(distribution), std::cend(distribution),
                                  std::cbegin(likelihoods), 0.0, std::plus<void> {},
                                  [k, n] (const auto p, const auto& haplotype_likelihoods) {
                                      return p * haplotype_likelihoods[k][n];
                                  });
    }
    
    template <std::size_t K>
    ResponsabilityVectors<K>
    init_responsabilities(const CompressedAlphas<K>& prior_alphas,
                          const ProbabilityVector& genotype_pobabilities,
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
                    ln_rho[k] = al[k] + expectation(genotype_pobabilities, read_likelihoods[s], k, n);
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
    template <std::size_t K>
    void update_responsabilities(ResponsabilityVectors<K>& result,
                                 const CompressedAlphas<K>& posterior_alphas,
                                 const ProbabilityVector& genotype_pobabilities,
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
                    ln_rho[k] = al[k] + expectation(genotype_pobabilities, read_likelihoods[s], k, n);
                }
                
                const auto ln_rho_norm = Maths::log_sum_exp(ln_rho);
                
                for (unsigned k {0}; k < K; ++k) {
                    result[s][n][k] = std::exp(ln_rho[k] - ln_rho_norm);
                }
            }
        }
    }
    
    template <std::size_t K>
    double sum(const ResponsabilityVector<K>& taus, const unsigned k)
    {
        return std::accumulate(std::cbegin(taus), std::cend(taus), 0.0,
                               [k] (const auto curr, const auto& tau) {
                                   return curr + tau[k];
                               });
    }
    
    template <std::size_t K>
    void update_alpha(CompressedAlpha<K>& alpha, const CompressedAlpha<K>& prior_alpha,
                      const ResponsabilityVector<K>& taus)
    {
        for (unsigned k {0}; k < K; ++k) {
            alpha[k] = prior_alpha[k] + sum(taus, k);
        }
    }
    
    template <std::size_t K>
    void update_alphas(CompressedAlphas<K>& alphas, const CompressedAlphas<K>& prior_alphas,
                       const ResponsabilityVectors<K>& responsabilities)
    {
        const auto S = alphas.size();
        assert(S == prior_alphas.size() && S == responsabilities.size());
        for (std::size_t s {0}; s < S; ++s) {
            update_alpha(alphas[s], prior_alphas[s], responsabilities[s]);
        }
    }
    
    template <std::size_t K>
    double marginalise(const ResponsabilityVectors<K>& responsabilities,
                       const CompressedReadLikelihoods<K>& read_likelihoods,
                       const std::size_t g)
    {
        double result {0};
        
        const auto S = read_likelihoods.size(); // num samples
        
        assert(S == responsabilities.size());
        
        for (std::size_t s {0}; s < S; ++s) {
            const auto N = read_likelihoods[s][0][0].size(); // num reads in sample s
            
            assert(responsabilities[s].size() == N);
            assert(responsabilities[s][0].size() == K);
            assert(read_likelihoods[s][g].size() == K);
            
            for (unsigned k {0}; k < K; ++k) {
                for (std::size_t n {0}; n < N; ++n) {
                    result += responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
                }
            }
        }
        
        return result;
    }
    
    template <std::size_t K>
    void update_genotype_log_posteriors(LogProbabilityVector& result,
                                        const LogProbabilityVector& genotype_log_priors,
                                        const ResponsabilityVectors<K>& responsabilities,
                                        const CompressedReadLikelihoods<K>& read_likelihoods)
    {
        for (std::size_t g {0}; g < result.size(); ++g) {
            result[g] = genotype_log_priors[g] + marginalise(responsabilities, read_likelihoods, g);
        }
        
        Maths::normalise_logs(result);
    }
    
    auto max_change(const CompressedAlpha<2>& lhs, const CompressedAlpha<2>& rhs)
    {
        return std::max(std::abs(lhs.front() - lhs.front()), std::abs(lhs.back() - lhs.back()));
    }
    
    auto max_change(const CompressedAlpha<3>& lhs, const CompressedAlpha<3>& rhs)
    {
        return std::max({std::abs(lhs[0] - rhs[0]), std::abs(lhs[1] - rhs[1]), std::abs(lhs[2] - rhs[2])});
    }
    
    template <std::size_t K>
    auto max_change(const CompressedAlpha<K>& lhs, const CompressedAlpha<K>& rhs)
    {
        double result {0};
        
        for (std::size_t k {0}; k < K; ++k) {
            const auto curr = std::abs(lhs[k] - rhs[k]);
            if (curr > result) result = curr;
        }
        
        return result;
    }
    
    template <std::size_t K>
    auto max_change(const CompressedAlphas<K>& prior_alphas,
                    const CompressedAlphas<K>& posterior_alphas)
    {
        assert(prior_alphas.size() == posterior_alphas.size());
        
        double result {0};
        
        for (std::size_t s {0}; s < prior_alphas.size(); ++s) {
            const auto curr = max_change(prior_alphas[s], posterior_alphas[s]);
            if (curr > result) result = curr;
        }
        
        return result;
    }
    
    template <std::size_t K>
    std::pair<bool, double> check_convergence(const CompressedAlphas<K>& prior_alphas,
                                              const CompressedAlphas<K>& posterior_alphas,
                                              const double prev_max_change,
                                              const double epsilon)
    {
        const auto new_max_change = max_change(prior_alphas, posterior_alphas);
        return std::make_pair(std::abs(new_max_change - prev_max_change) < epsilon, new_max_change);
    }
    
    // lower-bound calculation
    
    double expectation(const ProbabilityVector& genotype_posteriors,
                       const LogProbabilityVector& genotype_log_priors)
    {
        return std::inner_product(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                  std::cbegin(genotype_log_priors), 0.0);
    }
    
    template <std::size_t K>
    double dirichlet_expectation(const CompressedAlpha<K>& priors, const CompressedAlpha<K>& posteriors)
    {
        using boost::math::digamma;
        const auto da0 = digamma(sum(posteriors));
        return std::inner_product(std::cbegin(priors), std::cend(priors), std::cbegin(posteriors),
                                  0.0, std::plus<void> {},
                                  [da0] (const auto& prior, const auto& post) {
                                      return (prior - 1) * (digamma(post) - da0);
                                  }) - Maths::log_beta(priors);
    }
    
    template <std::size_t K>
    double expectation(const CompressedAlphas<K>& priors, const CompressedAlphas<K>& posteriors)
    {
        return std::inner_product(std::cbegin(priors), std::cend(priors), std::cbegin(posteriors),
                                  0.0, std::plus<void> {},
                                  [] (const auto& prior, const auto& post) {
                                      return dirichlet_expectation(prior, post);
                                  });
    }
    
    // E[ln p(Z_s | pi_s)]
    template <std::size_t K>
    double expectation(const ResponsabilityVector<K>& taus, const CompressedAlpha<K>& alpha)
    {
        using boost::math::digamma;
        
        const auto das = digamma(sum(alpha));
        
        double result {0};
        
        for (unsigned k {0}; k < K; ++k) {
            result += (digamma(alpha[k]) - das) * sum(taus, k);
        }
        
        return result;
    }
    
    // sum s E[ln p(Z_s | pi_s)]
    template <std::size_t K>
    double expectation(const ResponsabilityVectors<K>& taus, const CompressedAlphas<K>& alphas)
    {
        return std::inner_product(std::cbegin(taus), std::cend(taus), std::cbegin(alphas),
                                  0.0, std::plus<void> {},
                                  [] (const auto& tau, const auto& alpha) {
                                      return expectation(tau, alpha);
                                  });
    }
    
    template <std::size_t K>
    double expectation(const ResponsabilityVectors<K>& taus,
                       const CompressedReadLikelihoods<K>& log_likelihoods,
                       const std::size_t g)
    {
        double result {0};
        
        for (std::size_t s {0}; s < taus.size(); ++s) {
            for (std::size_t n {0}; n < taus[s].size(); ++n) {
                for (unsigned k {0}; k < K; ++k) {
                    result += taus[s][n][k] * log_likelihoods[s][g][k][n];
                }
            }
        }
        
        return result;
    }
    
    // E[ln p(R | Z, g)]
    template <std::size_t K>
    double expectation(const ProbabilityVector& genotype_posteriors,
                       const ResponsabilityVectors<K>& taus,
                       const CompressedReadLikelihoods<K>& log_likelihoods)
    {
        double result {0};
        
        for (std::size_t g {0}; g < genotype_posteriors.size(); ++g) {
            result += genotype_posteriors[g] * expectation(taus, log_likelihoods, g);
        }
        
        return result;
    }
    
    template <std::size_t K>
    double dirichlet_expectation(const CompressedAlpha<K>& posterior)
    {
        using boost::math::digamma;
        const auto da0 = digamma(sum(posterior));
        return std::accumulate(std::cbegin(posterior), std::cend(posterior), 0.0,
                               [da0] (const auto curr, const auto a) {
                                   return curr + ((a - 1) * (digamma(a) - da0));
                               }) - Maths::log_beta(posterior);
    }
    
    template <std::size_t K>
    double expectation(const CompressedAlphas<K>& posteriors)
    {
        return std::accumulate(std::cbegin(posteriors), std::cend(posteriors), 0.0,
                               [] (const auto curr, const auto& posterior) {
                                   return curr + dirichlet_expectation(posterior);
                               });
    }
    
    template <std::size_t K>
    double q_expectation(const Tau<K>& tau)
    {
        return std::accumulate(std::cbegin(tau), std::cend(tau), 0.0,
                               [] (const auto curr, const auto t) {
                                   return curr + (t * std::log(t));
                               });
    }
    
    // E [ln q(Z_s)]
    template <std::size_t K>
    double q_expectation(const ResponsabilityVector<K>& taus)
    {
        return std::accumulate(std::cbegin(taus), std::cend(taus), 0.0,
                               [] (const auto curr, const auto& tau) {
                                   return curr + q_expectation(tau);
                               });
    }
    
    // sum s E [ln q(Z_s)]
    template <std::size_t K>
    double q_expectation(const ResponsabilityVectors<K>& taus)
    {
        return std::accumulate(std::cbegin(taus), std::cend(taus), 0.0,
                               [] (const auto curr, const auto& t) {
                                   return curr + q_expectation(t);
                               });
    }
    
    template <std::size_t K>
    double calculate_lower_bound(const CompressedAlphas<K>& prior_alphas,
                                 const LogProbabilityVector& genotype_log_priors,
                                 const CompressedReadLikelihoods<K>& log_likelihoods,
                                 const CompressedLatents<K>& latents)
    {
        const auto& genotype_posteriors = latents.genotype_posteriors;
        const auto& genotype_log_posteriors = latents.genotype_log_posteriors;
        const auto& posterior_alphas = latents.alphas;
        const auto& taus = latents.responsabilities;
        
        double result {0};
        
        result += expectation(genotype_posteriors, genotype_log_priors);
        result += expectation(prior_alphas, posterior_alphas);
        result += expectation(taus, posterior_alphas);
        result += expectation(genotype_posteriors, taus, log_likelihoods);
        
        result -= expectation(genotype_posteriors, genotype_log_posteriors);
        result -= expectation(posterior_alphas);
        result -= q_expectation(taus);
        
        return result;
    }
    
    // Main algorithm - single seed
    
    // Starting iteration with given genotype_log_posteriors
    template <std::size_t K>
    CompressedLatents<K>
    run_variational_bayes(const CompressedAlphas<K>& prior_alphas,
                          const LogProbabilityVector& genotype_log_priors,
                          const CompressedReadLikelihoods<K>& log_likelihoods,
                          LogProbabilityVector genotype_log_posteriors,
                          const VariationalBayesParameters params)
    {
        assert(!prior_alphas.empty());
        assert(!genotype_log_priors.empty());
        assert(!log_likelihoods.empty());
        assert(prior_alphas.size() == log_likelihoods.size()); // num samples
        assert(log_likelihoods.front().size() == genotype_log_priors.size()); // num genotypes
        assert(params.max_iterations > 0);
        
        auto genotype_posteriors = exp(genotype_log_posteriors);
        
        auto posterior_alphas = prior_alphas;
        
        auto responsabilities = init_responsabilities<K>(posterior_alphas, genotype_posteriors,
                                                         log_likelihoods);
        
        assert(responsabilities.size() == log_likelihoods.size()); // num samples
        //assert(!responsabilities.front().empty());
        
        bool is_converged {false};
        double max_change {0};
        
        // main loop
        for (unsigned i {0}; i < params.max_iterations; ++i) {
            if (TRACE_MODE) {
                Logging::TraceLogger log {};
                stream(log) << "VB Iteration " << i;
            }
            
            update_genotype_log_posteriors(genotype_log_posteriors, genotype_log_priors,
                                           responsabilities, log_likelihoods);
            
            exp(genotype_posteriors, genotype_log_posteriors);
            
            update_alphas(posterior_alphas, prior_alphas, responsabilities);
            
            update_responsabilities(responsabilities, posterior_alphas, genotype_posteriors,
                                    log_likelihoods);
            
            std::tie(is_converged, max_change) = check_convergence(prior_alphas, posterior_alphas,
                                                                   max_change, params.epsilon);
            
            if (is_converged) {
                break;
            }
        }
        
        return {std::move(genotype_posteriors), std::move(genotype_log_posteriors),
                std::move(posterior_alphas), std::move(responsabilities)};
    }
    
    // Main algorithm - multiple seeds
    
    template <std::size_t K>
    std::pair<CompressedLatents<K>, double>
    run_variational_bayes(const CompressedAlphas<K>& prior_alphas,
                          const LogProbabilityVector& genotype_log_priors,
                          const CompressedReadLikelihoods<K>& log_likelihoods,
                          const VariationalBayesParameters params,
                          std::vector<LogProbabilityVector> seeds)
    {
        // Try the main algorithm from different seeds to check local optimum
        
        std::vector<CompressedLatents<K>> results;
        results.reserve(seeds.size());
        
        for (auto& seed : seeds) {
            results.emplace_back(run_variational_bayes(prior_alphas, genotype_log_priors,
                                                       log_likelihoods,
                                                       std::move(seed),
                                                       params));
        }
        
        std::vector<double> result_evidences(results.size());
        
        std::transform(std::cbegin(results), std::cend(results), std::begin(result_evidences),
                       [&] (const auto& latents) {
                           return calculate_lower_bound(prior_alphas, genotype_log_priors,
                                                        log_likelihoods, latents);
                       });
        
        const auto it = std::max_element(std::cbegin(result_evidences),
                                         std::cend(result_evidences));
        
        const auto idx = std::distance(std::cbegin(result_evidences), it);
        
        return std::make_pair(std::move(results[idx]), *it);
    }
    
    // Helpers
    
    CNV::Latents::GenotypeProbabilityMap
    expand(std::vector<Genotype<Haplotype>>&& genotypes, LogProbabilityVector&& genotype_log_posteriors)
    {
        CNV::Latents::GenotypeProbabilityMap result {};
        
        std::transform(std::make_move_iterator(std::begin(genotypes)),
                       std::make_move_iterator(std::end(genotypes)),
                       std::begin(genotype_log_posteriors),
                       std::inserter(result, std::begin(result)),
                       [] (auto&& g, auto p) {
                           return std::make_pair(std::move(g), p);
                       });
        
        return result;
    }
    
    template <std::size_t K>
    CNV::Latents::GenotypeMixturesDirichletAlphas expand(CompressedAlpha<K>& alpha)
    {
        return CNV::Latents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
    }
    
    template <std::size_t K>
    CNV::Latents::GenotypeMixturesDirichletAlphaMap
    expand(const std::vector<SampleIdType>& samples, CompressedAlphas<K>&& alphas)
    {
        CNV::Latents::GenotypeMixturesDirichletAlphaMap result {};
        
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                       std::inserter(result, std::begin(result)),
                       [] (const auto& sample, auto&& compressed_alpha) {
                           return std::make_pair(sample, expand(compressed_alpha));
                       });
        
        return result;
    }
    
    template <std::size_t K>
    CNV::InferredLatents
    expand(const std::vector<SampleIdType>& samples, std::vector<Genotype<Haplotype>>&& genotypes,
           CompressedLatents<K>&& inferred_latents, double evidence)
    {
        CNV::Latents posterior_latents {
            expand(std::move(genotypes), std::move(inferred_latents.genotype_posteriors)),
            expand(samples, std::move(inferred_latents.alphas))
        };
        
        return CNV::InferredLatents {std::move(posterior_latents), evidence};
    }
    
    LogProbabilityVector log_uniform_dist(const std::size_t n)
    {
        return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
    }
    
    auto generate_seeds(const std::vector<SampleIdType>& samples,
                        const std::vector<Genotype<Haplotype>>& genotypes,
                        const LogProbabilityVector& genotype_log_priors,
                        const CNV::Priors& priors,
                        const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
    {
        std::vector<LogProbabilityVector> result {};
        
        result.emplace_back(genotype_log_priors);
        result.emplace_back(log_uniform_dist(genotypes.size()));
        
        GenotypeModel::Individual germline_model {genotypes.front().ploidy(), priors.genotype_prior_model};
        
        for (const auto& sample : samples) {
            const auto latents = germline_model.infer_latents(sample, genotypes,
                                                              haplotype_log_likelihoods);
            result.emplace_back(latents.posteriors.genotype_probabilities);
            Maths::log_each(result.back());
        }
        
        return result;
    }
    
    // Main entry point
    
    template <std::size_t K>
    CNV::InferredLatents
    run_variational_bayes(const std::vector<SampleIdType>& samples,
                          std::vector<Genotype<Haplotype>>&& genotypes,
                          const CNV::Priors& priors,
                          const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                          const VariationalBayesParameters& params)
    {
        const auto prior_alphas = flatten_priors<K>(priors, samples);
        
        const auto genotype_log_priors = calculate_log_priors(genotypes, priors.genotype_prior_model);
        
        const auto log_likelihoods = compress<K>(genotypes, samples, haplotype_log_likelihoods);
        
        auto seeds = generate_seeds(samples, genotypes, genotype_log_priors,
                                    priors, haplotype_log_likelihoods);
        
        auto p = run_variational_bayes(prior_alphas, genotype_log_priors,
                                       log_likelihoods, params, std::move(seeds));
        
        return expand(samples, std::move(genotypes), std::move(p.first), p.second);
    }
    
    // Previously declared helpers
    
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
}
} // namespace GenotypeModel
} // namespace Octopus
