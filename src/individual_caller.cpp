//
//  individual_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "individual_caller.hpp"

#include <unordered_map>
#include <deque>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <utility>
#include <iostream>

#include <boost/iterator/zip_iterator.hpp>
#include <boost/tuple/tuple.hpp>

#include "genomic_region.hpp"
#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "maths.hpp"
#include "mappable_algorithms.hpp"
#include "read_utils.hpp"
#include "probability_matrix.hpp"
#include "coalescent_model.hpp"
#include "individual_genotype_model.hpp"
#include "germline_variant_call.hpp"
#include "reference_call.hpp"
#include "logging.hpp"

namespace Octopus
{
IndividualVariantCaller::CallerParameters::CallerParameters(double min_variant_posterior,
                                                            double min_refcall_posterior,
                                                            unsigned ploidy)
:
min_variant_posterior {min_variant_posterior},
min_refcall_posterior {min_refcall_posterior},
ploidy {ploidy}
{}

IndividualVariantCaller::IndividualVariantCaller(const ReferenceGenome& reference,
                                                 ReadPipe& read_pipe,
                                                 CandidateVariantGenerator&& candidate_generator,
                                                 VariantCaller::CallerParameters general_parameters,
                                                 CallerParameters specific_parameters)
:
VariantCaller {reference, read_pipe, std::move(candidate_generator), std::move(general_parameters)},
sample_ {read_pipe.get_samples().front()},
ploidy_ {specific_parameters.ploidy},
min_variant_posterior_ {specific_parameters.min_variant_posterior},
min_refcall_posterior_ {specific_parameters.min_refcall_posterior}
{}

// IndividualVariantCaller::Latents public methods

IndividualVariantCaller::Latents::Latents(const SampleIdType& sample,
                                          const std::vector<Haplotype>& haplotypes,
                                          std::vector<Genotype<Haplotype>>&& genotypes,
                                          ModelLatents&& latents)
:
genotype_posteriors_ {},
haplotype_posteriors_ {}
{
    GenotypeProbabilityMap genotype_posteriors {
        std::make_move_iterator(std::begin(genotypes)),
        std::make_move_iterator(std::end(genotypes))
    };
    
    insert_sample(sample, latents.genotype_probabilities, genotype_posteriors);
    
    genotype_posteriors_  = std::make_shared<GenotypeProbabilityMap>(std::move(genotype_posteriors));
    haplotype_posteriors_ = std::make_shared<HaplotypeProbabilityMap>(calculate_haplotype_posteriors(haplotypes));
}

std::shared_ptr<IndividualVariantCaller::Latents::HaplotypeProbabilityMap>
IndividualVariantCaller::Latents::get_haplotype_posteriors() const noexcept
{
    return haplotype_posteriors_;
}

std::shared_ptr<IndividualVariantCaller::Latents::GenotypeProbabilityMap>
IndividualVariantCaller::Latents::get_genotype_posteriors() const noexcept
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
    CoalescentModel prior_model {Haplotype {mapped_region(haplotypes.front()), reference_}};
    
    GenotypeModel::Individual model {ploidy_, prior_model};
    
    auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "There are " << genotypes.size() << " candidate genotypes";
    }
    
    const auto& sample = samples_.front();
    
    auto inferences = model.infer_latents(sample, genotypes, haplotype_likelihoods);
    
//    // TEST
//    GenotypeModel::Individual dummy_model {ploidy_ + 1, prior_model};
//    auto dummy_genotypes = generate_all_genotypes(haplotypes, ploidy_ + 1);
//    auto dummy_inferences = dummy_model.infer_latents(sample, dummy_genotypes, haplotype_likelihoods);
//    auto lp1 = std::log(0.99) + inferences.log_evidence;
//    auto lp2 = std::log(1.0 - 0.99) + dummy_inferences.log_evidence;
//    auto norm = Maths::log_sum_exp(lp1, lp2);
//    lp1 -= norm; lp2 -= norm;
//    if (lp1 < lp2) {
//        std::cout << inferences.log_evidence << " " << dummy_inferences.log_evidence << '\n';
//        auto it = std::max_element(std::cbegin(dummy_inferences.posteriors.genotype_probabilities),
//                                   std::cend(dummy_inferences.posteriors.genotype_probabilities));
//        auto itd = std::distance(std::cbegin(dummy_inferences.posteriors.genotype_probabilities), it);
//        ::debug::print_variant_alleles(dummy_genotypes[itd]);
//        std::cout << " " << *it << '\n';
//        std::cout << haplotypes.front().get_region() << std::endl;
//    }
//    // END TEST
    
    return std::make_unique<Latents>(sample, haplotypes, std::move(genotypes),
                                     std::move(inferences.posteriors));
}

namespace
{
using GM = GenotypeModel::Individual;

using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>::InnerMap;

using AllelePosteriorMap = std::unordered_map<Allele, double>;

struct VariantCall : public Mappable<VariantCall>
{
    VariantCall() = default;
    template <typename T>
    VariantCall(T&& variant, double posterior)
    :
    variant {std::forward<T>(variant)},
    posterior {posterior}
    {}
    
    GenomicRegion get_region() const { return variant.get_region(); }
    
    Variant variant;
    double posterior;
};

using VariantCalls = std::vector<VariantCall>;

struct GenotypeCall
{
    GenotypeCall() = default;
    template <typename T> GenotypeCall(T&& genotype, double posterior)
    :
    genotype {std::forward<T>(genotype)},
    posterior {posterior}
    {}
    
    Genotype<Allele> genotype;
    double posterior;
};

using GenotypeCalls = std::vector<GenotypeCall>;
} // namespace

namespace debug
{
    template <typename S>
    void print_genotype_posteriors(S&& stream, const GenotypeProbabilityMap& genotype_posteriors,
                                   std::size_t n = 5);
    void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                                   std::size_t n = 5);
    template <typename S>
    void print_allele_posteriors(S&& stream, const AllelePosteriorMap& allele_posteriors,
                                 std::size_t n = 10);
    void print_allele_posteriors(const AllelePosteriorMap& allele_posteriors,
                                 std::size_t n = 10);
//        void print_variant_calls(const VariantCallBlocks& calls);
} // namespace debug

namespace
{
// allele posterior calculations

using AlleleBools          = std::deque<bool>; // using std::deque because std::vector<bool> is evil
using GenotypeContainments = std::vector<AlleleBools>;

auto marginalise(const Allele& allele, const GenotypeProbabilityMap& genotype_posteriors)
{
    auto result = std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                                  0.0, [&allele] (const auto curr, const auto& p) {
                                      return curr + (contains(p.first, allele) ? p.second : 0);
                                  });
    
    if (result > 1.0) result = 1.0; // just to account for floating point error
    
    return result;
}

AllelePosteriorMap compute_allele_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                                             const std::vector<Allele>& alleles)
{
    AllelePosteriorMap result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        result.emplace(allele, marginalise(allele, genotype_posteriors));
    }
    
    return result;
}

// variant calling

VariantCalls call_variants(const std::vector<Variant>& candidates,
                           const AllelePosteriorMap& allele_posteriors,
                           const Genotype<Haplotype>& genotype_call,
                           const double min_posterior)
{
    assert(!allele_posteriors.empty());
    
    VariantCalls result {};
    result.reserve(candidates.size());
    
    for (const Variant& candidate : candidates) {
        const auto posterior = allele_posteriors.at(candidate.get_alt_allele());
        
        if (posterior >= min_posterior && contains_exact(genotype_call, candidate.get_alt_allele())) {
            result.emplace_back(candidate, posterior);
        }
    }
    
    result.shrink_to_fit();
    
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

double marginalise(const Genotype<Allele>& genotype, const GenotypeProbabilityMap& genotype_posteriors)
{
    return std::accumulate(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), 0.0,
                           [&genotype] (const double curr, const auto& p) {
                               return curr + ((contains(p.first, genotype)) ? p.second : 0.0);
                           });
}

GenotypeCalls call_genotypes(const Genotype<Haplotype>& genotype_call,
                             const GenotypeProbabilityMap& genotype_posteriors,
                             const std::vector<GenomicRegion>& variant_regions)
{
    GenotypeCalls result {};
    result.reserve(variant_regions.size());
    
    for (const auto& region : variant_regions) {
        auto spliced_genotype = splice<Allele>(genotype_call, region);
        
        const auto posterior = marginalise(spliced_genotype, genotype_posteriors);
        
        result.emplace_back(std::move(spliced_genotype), posterior);
    }
    
    return result;
}

// output

Octopus::VariantCall::GenotypeCall convert(GenotypeCall&& call)
{
    return Octopus::VariantCall::GenotypeCall {std::move(call.genotype), call.posterior};
}

std::unique_ptr<Octopus::VariantCall>
make_call(const SampleIdType& sample, VariantCall&& variant_call, GenotypeCall&& genotype_call)
{
    std::vector<std::pair<SampleIdType, Call::GenotypeCall>> tmp {
        std::make_pair(sample, convert(std::move(genotype_call)))
    };
    return std::make_unique<GermlineVariantCall>(std::move(variant_call.variant),
                                                 std::move(tmp), variant_call.posterior);
}

auto make_calls(const SampleIdType& sample, VariantCalls&& variant_calls,
                GenotypeCalls&& genotype_calls)
{
    std::vector<std::unique_ptr<Octopus::VariantCall>> result {};
    result.reserve(variant_calls.size());
    
    std::transform(std::make_move_iterator(std::begin(variant_calls)),
                   std::make_move_iterator(std::end(variant_calls)),
                   std::make_move_iterator(std::begin(genotype_calls)),
                   std::back_inserter(result),
                   [&sample] (VariantCall&& variant_call, GenotypeCall&& genotype_call) {
                       return make_call(sample, std::move(variant_call), std::move(genotype_call));
                   });
    
    return result;
}
} // namespace

std::vector<std::unique_ptr<Octopus::VariantCall>>
IndividualVariantCaller::call_variants(const std::vector<Variant>& candidates,
                                       const std::vector<Allele>& callable_alleles,
                                       CallerLatents* latents) const
{
    assert(!callable_alleles.empty());
    
    const auto dlatents = dynamic_cast<Latents*>(latents);
    
    const auto& genotype_posteriors = (*dlatents->genotype_posteriors_)[sample_];
    
    if (TRACE_MODE) {
        Logging::TraceLogger log {};
        debug::print_genotype_posteriors(stream(log), genotype_posteriors, -1);
    } else if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_genotype_posteriors(stream(log), genotype_posteriors);
    }
    
    const auto allele_posteriors = compute_allele_posteriors(genotype_posteriors, callable_alleles);
    
    if (TRACE_MODE) {
        Logging::TraceLogger log {};
        debug::print_allele_posteriors(stream(log), allele_posteriors, -1);
    } else if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        debug::print_allele_posteriors(stream(log), allele_posteriors);
    }
    
    const auto genotype_call = call_genotype(genotype_posteriors);
    
    auto variant_calls = Octopus::call_variants(candidates, allele_posteriors, genotype_call,
                                                min_variant_posterior_);
    
    const auto called_regions = extract_regions(variant_calls);
    
    auto genotype_calls = call_genotypes(genotype_call, genotype_posteriors, called_regions);
    
    return make_calls(sample_, std::move(variant_calls), std::move(genotype_calls));
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
        
        const GenomicRegion& get_region() const noexcept { return reference_allele.get_region(); }
        
        Allele reference_allele;
        double posterior;
    };
    
    using RefCalls = std::vector<RefCall>;
    
    double marginalise_reference_genotype(const Allele& reference_allele,
                                          const GenotypeProbabilityMap& sample_genotype_posteriors)
    {
        double result {0};
        
        for (const auto& genotype_posterior : sample_genotype_posteriors) {
            if (is_homozygous(genotype_posterior.first, reference_allele)) {
                result += genotype_posterior.second;
            }
        }
        
        return result;
    }
    
    RefCalls call_reference(const GenotypeProbabilityMap& genotype_posteriors,
                            const std::vector<Allele>& reference_alleles,
                            const ReadMap::mapped_type& reads, const double min_call_posterior)
    {
        RefCalls result {};
        
        if (reference_alleles.empty()) return result;
        
        result.reserve(reference_alleles.size());
        
        for (const auto& reference_allele : reference_alleles) {
            double posterior {0};
            
            if (has_coverage(reads, mapped_region(reference_allele))) {
                posterior = marginalise_reference_genotype(reference_allele,
                                                           genotype_posteriors);
            }
            
            if (posterior >= min_call_posterior) {
                result.emplace_back(reference_allele, posterior);
            }
        }
        
        result.shrink_to_fit();
        
        return result;
    }
} // namespace

std::vector<std::unique_ptr<Call>>
IndividualVariantCaller::call_reference(const std::vector<Allele>& alleles,
                                        CallerLatents* latents,
                                        const ReadMap& reads) const
{
    return {};
}

namespace debug
{
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
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    void print_genotype_posteriors(const GenotypeProbabilityMap& genotype_posteriors,
                                   const std::size_t n)
    {
        print_genotype_posteriors(std::cout, genotype_posteriors, n);
    }
    
    template <typename S>
    void print_allele_posteriors(S&& stream, const AllelePosteriorMap& allele_posteriors,
                                 const std::size_t n)
    {
        const auto m = std::min(n, allele_posteriors.size());
        
        if (m == allele_posteriors.size()) {
            stream << "Printing all allele posteriors " << '\n';
        } else {
            stream << "Printing top " << m << " allele posteriors " << '\n';
        }
        
        using AlleleReference = std::reference_wrapper<const Allele>;
        
        std::vector<std::pair<AlleleReference, double>> v {};
        v.reserve(allele_posteriors.size());
        
        std::copy(std::cbegin(allele_posteriors), std::cend(allele_posteriors),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&] (const auto& p) {
                          stream << p.first.get() << " " << p.second << '\n';
                      });
    }
    
    void print_allele_posteriors(const AllelePosteriorMap& allele_posteriors,
                                 const std::size_t n)
    {
        print_allele_posteriors(std::cout, allele_posteriors, n);
    }
} // namespace debug
} // namespace Octopus
