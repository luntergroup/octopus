//
//  individual_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "individual_caller.hpp"

#include <typeinfo>
#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <stdexcept>
#include <iostream>

#include "genomic_region.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "maths.hpp"
#include "mappable_algorithms.hpp"
#include "read_stats.hpp"
#include "probability_matrix.hpp"
#include "coalescent_model.hpp"
#include "individual.hpp"
#include "germline_variant_call.hpp"
#include "reference_call.hpp"
#include "logging.hpp"

#include "timers.hpp"

namespace octopus {

IndividualVariantCaller::IndividualVariantCaller(VariantCaller::Components&& components,
                                                 VariantCaller::Parameters general_parameters,
                                                 Parameters specific_parameters)
:
VariantCaller {std::move(components), std::move(general_parameters)},
parameters_ {std::move(specific_parameters)}
{
    if (parameters_.ploidy == 0) {
        throw std::logic_error {"IndividualVariantCaller: ploidy must be > 0"};
    }
}

IndividualVariantCaller::CallTypeSet IndividualVariantCaller::do_get_call_types() const
{
    return {std::type_index(typeid(GermlineVariantCall))};
}

// IndividualVariantCaller::Latents public methods

IndividualVariantCaller::Latents::Latents(const SampleName& sample,
                                          const std::vector<Haplotype>& haplotypes,
                                          std::vector<Genotype<Haplotype>>&& genotypes,
                                          ModelInferences&& inferences)
:
genotype_posteriors_ {},
haplotype_posteriors_ {},
model_log_evidence_ {inferences.log_evidence}
{
    GenotypeProbabilityMap genotype_posteriors {
        std::make_move_iterator(std::begin(genotypes)),
        std::make_move_iterator(std::end(genotypes))
    };
    
    insert_sample(sample, inferences.posteriors.genotype_probabilities, genotype_posteriors);
    
    genotype_posteriors_  = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes));
}

std::shared_ptr<IndividualVariantCaller::Latents::HaplotypeProbabilityMap>
IndividualVariantCaller::Latents::haplotype_posteriors() const noexcept
{
    return haplotype_posteriors_;
}

std::shared_ptr<IndividualVariantCaller::Latents::GenotypeProbabilityMap>
IndividualVariantCaller::Latents::genotype_posteriors() const noexcept
{
    return genotype_posteriors_;
}

// IndividualVariantCaller::Latents private methods

IndividualVariantCaller::Latents::HaplotypeProbabilityMap
IndividualVariantCaller::Latents::calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes)
{
    assert(genotype_posteriors_ != nullptr);
    
    HaplotypeProbabilityMap result {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        result.emplace(haplotype, 0.0);
    }
    
    const auto& sample = std::cbegin(*genotype_posteriors_)->first;
    
    for (const auto& p : (*genotype_posteriors_)[sample]) {
        for (const auto& haplotype : p.first.copy_unique_ref()) {
            result.at(haplotype) += p.second;
        }
    }
    
    return result;
}

std::unique_ptr<IndividualVariantCaller::CallerLatents>
IndividualVariantCaller::infer_latents(const std::vector<Haplotype>& haplotypes,
                                       const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    auto genotypes = generate_all_genotypes(haplotypes, parameters_.ploidy);
    
    if (debug_log_) stream(*debug_log_) << "There are " << genotypes.size() << " candidate genotypes";
    
    const auto prior_model = make_prior_model(haplotypes);
    
    const model::Individual model {prior_model, debug_log_};
    
    haplotype_likelihoods.prime(sample());
    
    auto inferences = model.infer_latents(genotypes, haplotype_likelihoods);
    
    return std::make_unique<Latents>(sample(), haplotypes, std::move(genotypes), std::move(inferences));
}

boost::optional<double>
IndividualVariantCaller::calculate_dummy_model_posterior(const std::vector<Haplotype>& haplotypes,
                                                         const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                                         const CallerLatents& latents) const
{
    return calculate_dummy_model_posterior(haplotypes, haplotype_likelihoods,
                                           dynamic_cast<const Latents&>(latents));
}

static auto calculate_dummy_model_posterior(const double normal_model_log_evidence,
                                            const double dummy_model_log_evidence)
{
    constexpr double normal_model_prior {0.9999999};
    constexpr double dummy_model_prior {1.0 - normal_model_prior};
    
    const auto normal_model_ljp = std::log(normal_model_prior) + normal_model_log_evidence;
    const auto dummy_model_ljp  = std::log(dummy_model_prior) + dummy_model_log_evidence;
    
    const auto norm = maths::log_sum_exp(normal_model_ljp, dummy_model_ljp);
    
    return std::exp(dummy_model_ljp - norm);
}

boost::optional<double>
IndividualVariantCaller::calculate_dummy_model_posterior(const std::vector<Haplotype>& haplotypes,
                                                         const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                                         const Latents& latents) const
{
    const auto genotypes = generate_all_genotypes(haplotypes, parameters_.ploidy + 1);
    
    const auto prior_model = make_prior_model(haplotypes);
    
    const model::Individual model {prior_model, debug_log_};
    
    haplotype_likelihoods.prime(sample());
    
    const auto inferences = model.infer_latents(genotypes, haplotype_likelihoods);
    
    return ::octopus::calculate_dummy_model_posterior(latents.model_log_evidence_,
                                                      inferences.log_evidence);
}

namespace
{
using GM = model::Individual;

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>::InnerMap;

using VariantReference  = std::reference_wrapper<const Variant>;
using VariantPosteriors = std::vector<std::pair<VariantReference, Phred<double>>>;

struct VariantCall : Mappable<VariantCall>
{
    VariantCall() = delete;
    VariantCall(const std::pair<VariantReference, Phred<double>>& p)
    : variant {p.first}, posterior {p.second} {}
    VariantCall(const Variant& variant, Phred<double> posterior)
    : variant {variant}, posterior {posterior} {}
    
    const GenomicRegion& mapped_region() const noexcept { return ::mapped_region(variant.get()); }
    
    VariantReference variant;
    Phred<double> posterior;
    bool is_dummy_filtered = false;
};

using VariantCalls = std::vector<VariantCall>;

struct GenotypeCall
{
    template <typename T> GenotypeCall(T&& genotype, Phred<double> posterior)
    : genotype {std::forward<T>(genotype)}, posterior {posterior} {}
    
    Genotype<Allele> genotype;
    Phred<double> posterior;
};

using GenotypeCalls = std::vector<GenotypeCall>;
} // namespace

namespace
{
// allele posterior calculations

auto compute_posterior(const Allele& allele, const GenotypeProbabilityMap& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                             0.0, [&allele] (const auto curr, const auto& p) {
                                 return curr + (contains(p.first, allele) ? 0.0 : p.second);
                             });
    return probability_to_phred(p);
}

VariantPosteriors compute_candidate_posteriors(const std::vector<Variant>& candidates,
                                               const GenotypeProbabilityMap& genotype_posteriors)
{
    VariantPosteriors result {};
    result.reserve(candidates.size());
    
    for (const auto& candidate : candidates) {
        result.emplace_back(candidate, compute_posterior(candidate.alt_allele(), genotype_posteriors));
    }
    
    return result;
}

// variant calling

bool contains_alt(const Genotype<Haplotype>& genotype_call, const VariantReference& candidate)
{
    return includes(genotype_call, candidate.get().alt_allele());
}

VariantCalls call_candidates(const VariantPosteriors& candidate_posteriors,
                             const Genotype<Haplotype>& genotype_call,
                             const Phred<double> min_posterior)
{
    VariantCalls result {};
    result.reserve(candidate_posteriors.size());
    
    std::copy_if(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
                 std::back_inserter(result),
                 [&genotype_call, min_posterior] (const auto& p) {
                     return p.second >= min_posterior && contains_alt(genotype_call, p.first);
                 });
    
    return result;
}

// variant genotype calling

auto call_genotype(const GenotypeProbabilityMap& genotype_posteriors)
{
    return std::max_element(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                            [] (const auto& lhs, const auto& rhs) {
                                return lhs.second < rhs.second;
                            })->first;
}

auto compute_posterior(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
{
    auto p = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                             [&genotype] (const double curr, const auto& p) {
                                 return curr + (contains(p.first, genotype) ? 0.0 : p.second);
                             });
    return probability_to_phred(p);
}

GenotypeCalls call_genotypes(const Genotype<Haplotype>& genotype_call,
                             const GenotypeProbabilityMap& genotype_posteriors,
                             const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    
    for (const auto& region : variant_regions) {
        auto spliced_genotype = splice<Allele>(genotype_call, region);
        
        const auto posterior = compute_posterior(spliced_genotype, genotype_posteriors);
        
        result.emplace_back(std::move(spliced_genotype), posterior);
    }
    
    return result;
}

// output

octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
{
    return octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<octopus::VariantCall>
transform_call(const SampleName& sample, VariantCall&& variant_call, GenotypeCall&& genotype_call)
{
    std::vector<std::pair<SampleName, Call::GenotypeCall>> tmp {
        std::make_pair(sample, convert(std::move(genotype_call)))
    };
    
    std::unique_ptr<octopus::VariantCall> result {
        std::make_unique<GermlineVariantCall>(variant_call.variant.get(), std::move(tmp),
                                              variant_call.posterior)
    };
    
    return result;
}

auto transform_calls(const SampleName& sample, VariantCalls&& variant_calls,
                     GenotypeCalls&& genotype_calls)
{
    std::vector<std::unique_ptr<octopus::VariantCall>> result {};
    result.reserve(variant_calls.size());
    
    std::transform(std::make_move_iterator(std::begin(variant_calls)),
                   std::make_move_iterator(std::end(variant_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)),
                   std::back_inserter(result),
                   [&sample] (VariantCall&& variant_call, GenotypeCall&& genotype_call) {
                       return transform_call(sample, std::move(variant_call), std::move(genotype_call));
                   });
    
    return result;
}
} // namespace

std::vector<std::unique_ptr<octopus::VariantCall>>
IndividualVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                       const CallerLatents& latents) const
{
    return call_variants(candidates, dynamic_cast<const Latents&>(latents));
}
    
namespace debug {
    void log(const GenotypeProbabilityMap& genotype_posteriors,
             boost::optional<logging::DebugLogger>& debug_log,
             boost::optional<logging::TraceLogger>& trace_log);
    void log(const VariantPosteriors& candidate_posteriors,
             boost::optional<logging::DebugLogger>& debug_log,
             boost::optional<logging::TraceLogger>& trace_log);
    void log(const Genotype<Haplotype>& called_genotype,
             boost::optional<logging::DebugLogger>& debug_log);
} // namespace debug

std::vector<std::unique_ptr<octopus::VariantCall>>
IndividualVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                       const Latents& latents) const
{
    const auto& genotype_posteriors = (*latents.genotype_posteriors_)[sample()];
    
    debug::log(genotype_posteriors, debug_log_, trace_log_);
    
    const auto candidate_posteriors = compute_candidate_posteriors(candidates, genotype_posteriors);
    
    debug::log(candidate_posteriors, debug_log_, trace_log_);
    
    const auto genotype_call = call_genotype(genotype_posteriors);
    
    auto variant_calls = call_candidates(candidate_posteriors, genotype_call,
                                         parameters_.min_variant_posterior);
    
    const auto called_regions = extract_regions(variant_calls);
    
    auto genotype_calls = call_genotypes(genotype_call, genotype_posteriors, called_regions);
    
    return transform_calls(sample(), std::move(variant_calls), std::move(genotype_calls));
}

namespace
{
    // reference genotype calling
    
    struct RefCall : public Mappable<RefCall>
    {
        RefCall() = default;
        
        template <typename A>
        RefCall(A&& reference_allele, double posterior)
        :
        reference_allele {std::forward<A>(reference_allele)},
        posterior {posterior}
        {}
        
        const GenomicRegion& mapped_region() const noexcept { return reference_allele.mapped_region(); }
        
        Allele reference_allele;
        double posterior;
    };
    
    using RefCalls = std::vector<RefCall>;
    
//    double marginalise_reference_genotype(const Allele& reference_allele,
//                                          const GenotypeProbabilityMap& sample_genotype_posteriors)
//    {
//        double result {0};
//        
//        for (const auto& genotype_posterior : sample_genotype_posteriors) {
//            if (is_homozygous(genotype_posterior.first, reference_allele)) {
//                result += genotype_posterior.second;
//            }
//        }
//        
//        return result;
//    }
    
//    RefCalls call_reference(const GenotypeProbabilityMap& genotype_posteriors,
//                            const std::vector<Allele>& reference_alleles,
//                            const ReadMap::mapped_type& reads, const double min_call_posterior)
//    {
//        RefCalls result {};
//        
//        if (reference_alleles.empty()) return result;
//        
//        result.reserve(reference_alleles.size());
//        
//        for (const auto& reference_allele : reference_alleles) {
//            double posterior {0};
//            
//            if (has_coverage(reads, mapped_region(reference_allele))) {
//                posterior = marginalise_reference_genotype(reference_allele,
//                                                           genotype_posteriors);
//            }
//            
//            if (posterior >= min_call_posterior) {
//                result.emplace_back(reference_allele, posterior);
//            }
//        }
//        
//        result.shrink_to_fit();
//        
//        return result;
//    }
} // namespace

std::vector<std::unique_ptr<ReferenceCall>>
IndividualVariantCaller::call_reference(const std::vector<Allele>& alleles,
                                        const CallerLatents& latents,
                                        const ReadMap& reads) const
{
    return call_reference(alleles, dynamic_cast<const Latents&>(latents), reads);
}

std::vector<std::unique_ptr<ReferenceCall>>
IndividualVariantCaller::call_reference(const std::vector<Allele>& alleles, const Latents& latents,
                                        const ReadMap& reads) const
{
    return {};
}

const SampleName& IndividualVariantCaller::sample() const noexcept
{
    return samples_.front();
}

CoalescentModel IndividualVariantCaller::make_prior_model(const std::vector<Haplotype>& haplotypes) const
{
    return CoalescentModel {
        Haplotype {mapped_region(haplotypes.front()), reference_},
        parameters_.snp_heterozygosity,
        parameters_.indel_heterozygosity
    };
}

namespace debug {
    template <typename S>
    void print_genotype_posteriors(S&& stream,
                                   const GenotypeProbabilityMap& genotype_posteriors,
                                   const std::size_t n)
    {
        const auto m = std::min(n, genotype_posteriors.size());
        
        if (m == genotype_posteriors.size()) {
            stream << "Printing all genotype posteriors " << '\n';
        } else {
            stream << "Printing top " << m << " genotype posteriors " << '\n';
        }
        
        using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
        
        std::vector<std::pair<GenotypeReference, double>> v {};
        v.reserve(genotype_posteriors.size());
        
        std::copy(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                                   const std::size_t n)
    {
        print_genotype_posteriors(std::cout, genotype_posteriors, n);
    }
    
    template <typename S>
    void print_candidate_posteriors(S&& stream, const VariantPosteriors& candidate_posteriors,
                                    const std::size_t n)
    {
        const auto m = std::min(n, candidate_posteriors.size());
        
        if (m == candidate_posteriors.size()) {
            stream << "Printing all candidate variant posteriors " << '\n';
        } else {
            stream << "Printing top " << m << " candidate variant posteriors " << '\n';
        }
        
        std::vector<std::pair<VariantReference, Phred<double>>> v {};
        v.reserve(candidate_posteriors.size());
        
        std::copy(std::cbegin(candidate_posteriors), std::cend(candidate_posteriors),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          stream << p.first.get() << " " << p.second.probability_true() << '\n';
                      });
    }
    
    void print_candidate_posteriors(const VariantPosteriors& candidate_posteriors,
                                    const std::size_t n)
    {
        print_candidate_posteriors(std::cout, candidate_posteriors, n);
    }
    
    void log(const GenotypeProbabilityMap& genotype_posteriors,
             boost::optional<logging::DebugLogger>& debug_log,
             boost::optional<logging::TraceLogger>& trace_log)
    {
        if (trace_log) {
            print_genotype_posteriors(stream(*trace_log), genotype_posteriors, -1);
        }
        if (debug_log) {
            print_genotype_posteriors(stream(*debug_log), genotype_posteriors, 5);
        }
    }
    
    void log(const VariantPosteriors& candidate_posteriors,
             boost::optional<logging::DebugLogger>& debug_log,
             boost::optional<logging::TraceLogger>& trace_log)
    {
        if (trace_log) {
            print_candidate_posteriors(stream(*trace_log), candidate_posteriors, -1);
        }
        if (debug_log) {
            print_candidate_posteriors(stream(*debug_log), candidate_posteriors, 5);
        }
    }
} // namespace debug
} // namespace octopus
