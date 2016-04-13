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

#include "fixed_ploidy_genotype_likelihood_model.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
    // public methods
    
    Cancer::Cancer(std::vector<SampleIdType> samples, const unsigned ploidy, Priors priors)
    :
    Cancer {std::move(samples), ploidy, std::move(priors), AlgorithmParameters {}}
    {}
        
    Cancer::InferredLatents::InferredLatents(Latents&& posterior_latents, double approx_log_evidence)
    :
    posterior_latents {std::move(posterior_latents)},
    approx_log_evidence {approx_log_evidence}
    {}
    
    Cancer::Cancer(std::vector<SampleIdType> samples, const unsigned ploidy, Priors priors,
                   AlgorithmParameters parameters)
    :
    samples_ {std::move(samples)},
    ploidy_ {ploidy},
    priors_ {std::move(priors)},
    parameters_ {parameters}
    {}
        
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
    Cancer::InferredLatents
    run_variational_bayes(const std::vector<SampleIdType>& samples,
                          std::vector<CancerGenotype<Haplotype>>&& genotypes,
                          const Cancer::Priors& priors,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods,
                          const VariationalBayesParameters& params);
    
    // Cancer public
    
    Cancer::InferredLatents
    Cancer::infer_latents(std::vector<CancerGenotype<Haplotype>> genotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        assert(!genotypes.empty());
        
        if (DEBUG_MODE) {
            Logging::DebugLogger log {};
            stream(log) << "There are " << genotypes.size() << " initial cancer genotypes";
        }
        
        const VariationalBayesParameters vb_params {parameters_.epsilon, parameters_.max_iterations};
        
        assert(ploidy_ < 3);
        
        if (ploidy_ == 1) {
            return run_variational_bayes<2>(samples_, std::move(genotypes), priors_,
                                            haplotype_likelihoods, vb_params);
        }
        
        return run_variational_bayes<3>(samples_, std::move(genotypes), priors_,
                                        haplotype_likelihoods, vb_params);
    }
    
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
    struct CompressedLatents
    {
        ProbabilityVector genotype_posteriors;
        CompressedAlphas<K> alphas;
    };
    
    template <std::size_t K>
    using Tau = std::array<double, K>;
    
    template <std::size_t K>
    using ResponsabilityVector = std::vector<Tau<K>>;
    
    template <std::size_t K>
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
    
    auto bundle(const Haplotype& a, const Haplotype& b)
    {
        return std::array<std::reference_wrapper<const Haplotype>, 2> {std::cref(a), std::cref(b)};
    }
    
    template <std::size_t K>
    CompressedAlpha<K> compress(const Cancer::Priors::GenotypeMixturesDirichletAlphas& alpha)
    {
        CompressedAlpha<K> result;
        std::copy_n(std::cbegin(alpha), K, std::begin(result));
        return result;
    }
    
    template <std::size_t K>
    CompressedAlphas<K> flatten_priors(const Cancer::Priors& priors, const std::vector<SampleIdType>& samples)
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
    
    template <std::size_t K>
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
    
    template <std::size_t K>
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
                                  std::cbegin(likelihoods), 0.0,
                                  std::plus<void> {},
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
                double curr {0}; // DEBUG
                for (std::size_t n {0}; n < N; ++n) {
//                    if (TRACE_MODE) {
//                        Logging::TraceLogger log {};
//                        if (g == 82 || g == 90) {
//                            stream(log) << "s: " << s << " n: " << n << " k: " << k << " : "
//                                    << responsabilities[s][n][k] << " * " << read_likelihoods[s][g][k][n]
//                                    << " = "<< responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
//                        }
//                    }
                    curr += responsabilities[s][n][k] * read_likelihoods[s][g][k][n]; // DEBUG
                    result += responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
                }
//                if (TRACE_MODE) {
//                    Logging::TraceLogger log {};
//                    if (g == 82 || g == 90) {
//                        stream(log) << k << " total = " << curr;
//                    }
//                }
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
//            if (TRACE_MODE) {
//                Logging::TraceLogger log {};
//                if (g == 82 || g == 90) {
//                    stream(log) << "g = " << g;
//                }
//            }
            result[g] = genotype_log_priors[g] + marginalise(responsabilities, read_likelihoods, g);
//            if (TRACE_MODE) {
//                Logging::TraceLogger log {};
//                if (g == 82 || g == 90) {
//                    stream(log) << "result " << result[g];
//                }
//            }
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
        assert(!responsabilities.front().empty());
        
        bool is_converged {false};
        double max_change {0};
        
        // main loop
        for (unsigned i {0}; i < params.max_iterations; ++i) {
            if (TRACE_MODE) {
                Logging::TraceLogger log {};
                stream(log) << "VB Iteration " << i;
            }
            //std::cout << "VB Iteration " << i << '\n';
            
            update_genotype_log_posteriors(genotype_log_posteriors, genotype_log_priors,
                                           responsabilities, log_likelihoods);
            
            exp(genotype_posteriors, genotype_log_posteriors);
            
            update_alphas(posterior_alphas, prior_alphas, responsabilities);
            
            update_responsabilities(responsabilities, posterior_alphas, genotype_posteriors,
                                    log_likelihoods);
            
            std::tie(is_converged, max_change) = check_convergence(prior_alphas, posterior_alphas,
                                                                   max_change, params.epsilon);
            
            //std::cout << "max change = " << max_change << std::endl;
            if (is_converged) {
                break;
            }
        }
        
        return {std::move(genotype_posteriors), std::move(posterior_alphas)};
    }
    
    // Helpers
    
    LogProbabilityVector log_uniform_dist(const std::size_t n)
    {
        return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
    }
    
    // Evidence calculation
    
    template <std::size_t K>
    auto expectation(const CompressedAlpha<K>& alpha)
    {
        const auto a0 = sum(alpha);
        std::array<double, K> result;
        std::transform(std::cbegin(alpha), std::cend(alpha), std::begin(result),
                       [a0] (const auto a) { return a / a0; });
        return result;
    }
    
    template <std::size_t K>
    auto log_marginal(const CompressedGenotype<K>& log_likelihoods,
                      const CompressedAlpha<K>& alpha,
                      std::array<double, K> pi,
                      const std::size_t n)
    {
        Maths::log_each(pi);
        
        for (unsigned k {0}; k < K; ++k) {
            pi[k] += log_likelihoods[k][n];
        }
        
        return Maths::log_sum_exp(pi);
    }
    
    template <std::size_t K>
    auto log_marginal(const CompressedGenotype<K>& log_likelihoods,
                      const CompressedAlpha<K>& alpha, std::array<double, K> pi)
    {
        auto result = Maths::log_dirichlet(alpha, pi);
        
        const auto N = log_likelihoods[0].size();
        
        for (std::size_t n {0}; n < N; ++n) {
            result += log_marginal(log_likelihoods, alpha, pi, n);
        }
        
        return result;
    }
    
    template <std::size_t K>
    double approx_log_marginal(const CompressedGenotype<K>& log_likelihoods,
                               const CompressedAlpha<K>& alpha)
    {
        const auto pi = expectation(alpha);
        return log_marginal(log_likelihoods, alpha, pi);
    }
    
    template <std::size_t K>
    double log_marginal(const CompressedReadLikelihoods<K>& log_likelihoods,
                        const CompressedAlphas<K>& alphas,
                        const std::size_t g)
    {
        double result {0};
        
        const auto S = alphas.size();
        
        for (std::size_t s {0}; s < S; ++s) {
            result += approx_log_marginal(log_likelihoods[s][g], alphas[s]);
        }
        
        return result;
    }
    
    template <std::size_t K>
    auto approx_log_evidence(const CompressedReadLikelihoods<K>& log_likelihoods,
                             const CompressedLatents<K>& latents)
    {
        const auto& genotype_posteriors = latents.genotype_posteriors;
        const auto& alphas = latents.alphas;
        
        const auto G = genotype_posteriors.size();
        
        std::vector<double> log_jcs(G);
        
        for (std::size_t g {0}; g < G; ++g) {
            log_jcs[g] = std::log(genotype_posteriors[g]) + log_marginal(log_likelihoods, alphas, g);
        }
        
        return Maths::log_sum_exp(log_jcs);
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
                       [&log_likelihoods] (const auto& latents) {
                           return approx_log_evidence(log_likelihoods, latents);
                       });
        
        const auto it = std::max_element(std::cbegin(result_evidences),
                                         std::cend(result_evidences));
        
        const auto idx = std::distance(std::cbegin(result_evidences), it);
        
        return std::make_pair(std::move(results[idx]), *it);
    }
    
    // Helpers
    
    Cancer::Latents::GenotypeProbabilityMap
    expand(std::vector<CancerGenotype<Haplotype>>&& genotypes,
           LogProbabilityVector&& genotype_log_posteriors)
    {
        Cancer::Latents::GenotypeProbabilityMap result {};
        
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
    Cancer::Latents::GenotypeMixturesDirichletAlphas expand(CompressedAlpha<K>& alpha)
    {
        return Cancer::Latents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
    }
    
    template <std::size_t K>
    Cancer::Latents::GenotypeMixturesDirichletAlphaMap
    expand(const std::vector<SampleIdType>& samples, CompressedAlphas<K>&& alphas)
    {
        Cancer::Latents::GenotypeMixturesDirichletAlphaMap result {};
        
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                       std::inserter(result, std::begin(result)),
                       [] (const auto& sample, auto&& compressed_alpha) {
                           return std::make_pair(sample, expand(compressed_alpha));
                       });
        
        return result;
    }
    
    template <std::size_t K>
    Cancer::InferredLatents
    expand(const std::vector<SampleIdType>& samples, std::vector<CancerGenotype<Haplotype>>&& genotypes,
           CompressedLatents<K>&& inferred_latents, const double evidence)
    {
        Cancer::Latents posterior_latents {
            expand(std::move(genotypes), std::move(inferred_latents.genotype_posteriors)),
            expand(samples, std::move(inferred_latents.alphas))
        };
        
        return Cancer::InferredLatents {std::move(posterior_latents), evidence};
    }
        
    auto infer_with_germline_model(const SampleIdType& sample,
                                   const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                   const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                                   const SomaticModel& genotype_prior_model)
    {
        assert(!genotypes.empty());
        
        const auto ploidy = genotypes.front().ploidy();
        
        FixedPloidyGenotypeLikelihoodModel likelihood_model {ploidy, haplotype_log_likelihoods};
        
        std::vector<double> result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&] (const auto& genotype) {
                           return std::log(genotype_prior_model.evaluate(genotype))
                                + likelihood_model.log_likelihood(sample, genotype.get_germline_genotype());
                       });
        
        Maths::normalise_logs(result);
        
        return result;
    }
    
    auto generate_seeds(const std::vector<SampleIdType>& samples,
                        const std::vector<CancerGenotype<Haplotype>>& genotypes,
                        const LogProbabilityVector& genotype_log_priors,
                        const Cancer::Priors& priors,
                        const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
    {
        std::vector<LogProbabilityVector> result {};
        result.reserve(2 + samples.size());
        
        result.emplace_back(genotype_log_priors);
        result.emplace_back(log_uniform_dist(genotypes.size()));
        
        for (const auto& sample : samples) {
            result.emplace_back(infer_with_germline_model(sample, genotypes,
                                                          haplotype_log_likelihoods,
                                                          priors.genotype_prior_model));
        }
        
        return result;
    }
    
    // Main entry point
    
    template <std::size_t K>
    Cancer::InferredLatents
    run_variational_bayes(const std::vector<SampleIdType>& samples,
                          std::vector<CancerGenotype<Haplotype>>&& genotypes,
                          const Cancer::Priors& priors,
                          const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                          const VariationalBayesParameters& params)
    {
        const auto prior_alphas = flatten_priors<K>(priors, samples);
        
        const auto genotype_log_priors = calculate_log_priors(genotypes, priors.genotype_prior_model);
        
//        if (TRACE_MODE) {
//            Logging::TraceLogger log {};
//            auto slog = stream(log);
//            slog << "All candidate cancer genotypes: " << '\n';
//            for (std::size_t g {0}; g < genotypes.size(); ++g) {
//                slog << g << " ";
//                debug::print_variant_alleles(slog, genotypes[g]);
//                slog << " " << genotype_log_priors[g] << '\n';
//            }
//        }
        
        const auto log_likelihoods = compress<K>(genotypes, samples, haplotype_log_likelihoods);
        
        auto seeds = generate_seeds(samples, genotypes, genotype_log_priors,
                                    priors, haplotype_log_likelihoods);
        
        auto p = run_variational_bayes(prior_alphas, genotype_log_priors,
                                       log_likelihoods, params, seeds);
        
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
    
    namespace debug
    {
        
    } // namespace debug
    
    } // namespace GenotypeModel
} // namespace Octopus
