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
    
    Cancer::InferredLatents
    Cancer::infer_latents(const std::vector<Haplotype>& haplotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        auto genotypes = generate_all_cancer_genotypes(haplotypes, ploidy_);
        
        return infer_latents(std::move(genotypes), haplotype_likelihoods);
    }
    
    // Compressed types used by the Variational Bayes model
    
    template <std::size_t K>
    using CompressedAlpha = std::array<double, K>;
    
    template <std::size_t K>
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
    
    template <std::size_t K>
    using CompressedGenotype = std::array<ReadLikelihoods, K>;
    template <std::size_t K>
    using CompressedGenotypes = std::vector<CompressedGenotype<K>>;
    template <std::size_t K>
    using CompressedReadLikelihoods = std::vector<CompressedGenotypes<K>>;
    
    template <std::size_t K>
    struct CompressedInferredLatents
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
    
    template <typename Container>
    void normalise_logs(Container& logs)
    {
        const auto norm = Maths::log_sum_exp(logs);
        for (auto& p : logs) p -= norm;
    }
    
    auto bundle(const Haplotype& a, const Haplotype& b)
    {
        return std::array<std::reference_wrapper<const Haplotype>, 2> {std::cref(a), std::cref(b)};
    }
    
    double probability_of_somatic(const Haplotype& somatic_haplotype,
                                  const Haplotype& germline_haplotype)
    {
        CoalescentModel model {germline_haplotype};
        return model.evaluate(bundle(germline_haplotype, somatic_haplotype));
    }
    
    // p(somatic | germline) = 1 / M sum k = 1 -> M p(somatic | germline_k) (M = germline ploidy)
    double probability_of_somatic(const Haplotype& somatic_haplotype,
                                  const Genotype<Haplotype>& germline_genotype)
    {
        const auto norm = 1.0 / germline_genotype.ploidy();
        
        return std::accumulate(std::cbegin(germline_genotype), std::cend(germline_genotype),
                               0.0,
                               [&somatic_haplotype] (const auto curr,
                                                     const Haplotype& germline_haplotype) {
                                   return curr + probability_of_somatic(somatic_haplotype,
                                                                        germline_haplotype);
                               }) / norm;
    }

    double calculate_log_prior(const CancerGenotype<Haplotype>& genotype,
                               const CoalescentModel& germline_prior_model)
    {
        const auto& germline = genotype.get_germline_genotype();
        const auto& somatic  = genotype.get_cancer_element();
        
        const auto germline_prior = germline_prior_model.evaluate(germline);
        
        const auto somatic_probability_given_germline = probability_of_somatic(somatic, germline);
        
        return std::log(germline_prior) + std::log(somatic_probability_given_germline);
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
                    ln_rho[k] = al[k] + expectation(genotype_priors, read_likelihoods[s], k, n);
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
                    ln_rho[k] = al[k] + expectation(genotype_posteriors, read_likelihoods[s], k, n);
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
                double curr {0};
                
                for (std::size_t n {0}; n < N; ++n) {
                    if (TRACE_MODE) {
                        Logging::TraceLogger log {};
                        if (g == 125 || g == 210) {
                            stream(log) << "n: " << n << " k: " << k << " : "
                                    << responsabilities[s][n][k] << " * " << read_likelihoods[s][g][k][n]
                                    << " = "<< responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
                        }
                    }
                    curr += responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
                    result += responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
                }
                
                if (TRACE_MODE) {
                    Logging::TraceLogger log {};
                    if (g == 125 || g == 210) {
                        stream(log) << k << " total = " << curr;
                    }
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
            if (TRACE_MODE) {
                Logging::TraceLogger log {};
                if (g == 125 || g == 210) {
                    stream(log) << "g = " << g;
                }
            }
            
            result[g] = genotype_log_priors[g] + marginalise(responsabilities, read_likelihoods, g);
            
            if (TRACE_MODE) {
                Logging::TraceLogger log {};
                if (g == 125 || g == 210) {
                    stream(log) << "result " << result[g];
                }
            }
        }
        normalise_logs(result);
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
    
    // The main algorithm
    template <std::size_t K>
    CompressedInferredLatents<K>
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
    
    LogProbabilityVector log_uniform_dist(const std::size_t n)
    {
        return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
    }
    
    template <std::size_t K>
    auto log_likelihood(const std::size_t g,
                        const std::vector<std::array<double, K>>& pis,
                        const CompressedAlphas<K>& alphas,
                        const LogProbabilityVector& genotype_log_probabilities,
                        const CompressedReadLikelihoods<K>& log_likelihoods)
    {
        auto result = genotype_log_probabilities[g];
        
        for (std::size_t s {0}; s < log_likelihoods.size(); ++s) {
            result += Maths::log_dirichlet(alphas[s], pis[s]);
            
            const auto N = log_likelihoods[s][0][0].size();
            for (std::size_t n {0}; n < N; ++n) {
                double tmp {0};
                for (unsigned k {0}; k < K; ++k) {
                    tmp += pis[s][k] * std::exp(log_likelihoods[s][g][k][n]);
                }
                result += std::log(tmp);
            }
        }
        
        return result;
    }
    
    template <std::size_t K>
    CompressedInferredLatents<K>
    run_variational_bayes(const CompressedAlphas<K>& prior_alphas,
                          const LogProbabilityVector& genotype_log_priors,
                          const CompressedReadLikelihoods<K>& log_likelihoods,
                          const VariationalBayesParameters params)
    {
        // Try the main algorithm from different seeds to check local optimum
        
        auto result1 = run_variational_bayes(prior_alphas, genotype_log_priors,
                                            log_likelihoods, genotype_log_priors,
                                            params);
        
        auto result2 = run_variational_bayes(prior_alphas, genotype_log_priors,
                                             log_likelihoods,
                                             log_uniform_dist(genotype_log_priors.size()),
                                             params);
        
        return result1;
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
    
    template <std::size_t K>
    Cancer::InferredLatents::GenotypeMixturesDirichletAlphas expand(CompressedAlpha<K>& alpha)
    {
        return Cancer::InferredLatents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
    }
    
    template <std::size_t K>
    Cancer::InferredLatents::GenotypeMixturesDirichletAlphaMap
    expand(const std::vector<SampleIdType>& samples, CompressedAlphas<K>&& alphas)
    {
        Cancer::InferredLatents::GenotypeMixturesDirichletAlphaMap result {};
        
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                       std::inserter(result, std::begin(result)),
                       [] (const auto& sample, auto&& compressed_alpha) {
                           return std::make_pair(sample, expand(compressed_alpha));
                       });
        
        return result;
    }
    
    template <std::size_t K>
    Cancer::InferredLatents
    expand(const std::vector<SampleIdType>& samples,
           std::vector<CancerGenotype<Haplotype>>&& genotypes,
           CompressedInferredLatents<K>&& inferred_latents)
    {
        return Cancer::InferredLatents {
            expand(std::move(genotypes), std::move(inferred_latents.genotype_posteriors)),
            expand(samples, std::move(inferred_latents.alphas))
        };
    }
    
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
        
        if (TRACE_MODE) {
            Logging::TraceLogger log {};
            auto slog = stream(log);
            slog << "All candidate cancer genotypes: " << '\n';
            for (std::size_t g {0}; g < genotypes.size(); ++g) {
                slog << g << " ";
                debug::print_variant_alleles(slog, genotypes[g]);
                slog << " " << genotype_log_priors[g] << '\n';
            }
        }
        
        const auto log_likelihoods = compress<K>(genotypes, samples, haplotype_log_likelihoods);
        
        auto inferred_latents = run_variational_bayes(prior_alphas,
                                                      genotype_log_priors,
                                                      log_likelihoods,
                                                      params);
        
        return expand(samples, std::move(genotypes), std::move(inferred_latents));
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
