//
//  variant_call.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef call_hpp
#define call_hpp

#include <vector>
#include <unordered_map>
#include <utility>
#include <iterator>

#include <boost/optional.hpp>

#include "common.hpp"
#include "phred.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
#include "allele.hpp"
#include "genotype.hpp"
#include "vcf_record.hpp"
#include "reference_genome.hpp"

namespace Octopus
{
    class Call : public Mappable<Call>
    {
    public:
        class PhaseCall
        {
        public:
            PhaseCall() = delete;
            
            template <typename R> PhaseCall(R&& region, Phred<double> score);
            
            const GenomicRegion& region() const noexcept { return region_; }
            
            Phred<double> score() const noexcept { return score_; };
            
        private:
            GenomicRegion region_;
            Phred<double> score_;
        };
        
        struct GenotypeCall
        {
            GenotypeCall() = delete;
            
            template <typename G> GenotypeCall(G&& genotype, Phred<double> posterior);
            template <typename G, typename P> GenotypeCall(G&& genotype, Phred<double> posterior, P&& phase);
            
            Genotype<Allele> genotype;
            Phred<double> posterior;
            boost::optional<PhaseCall> phase;
        };
        
        Call() = delete;
        
        explicit Call(Phred<double> quality);
        
        template <typename T> explicit Call(T&& genotype_calls, Phred<double> quality);
        
        virtual ~Call() = default;
        
        Call(const Call&)            = default;
        Call& operator=(const Call&) = default;
        Call(Call&&)                 = default;
        Call& operator=(Call&&)      = default;
        
        Phred<double> quality() const noexcept;
        
        GenotypeCall& get_genotype_call(const SampleIdType& sample);
        const GenotypeCall& get_genotype_call(const SampleIdType& sample) const;
        
        bool is_phased(const SampleIdType& sample) const;
        bool all_phased() const noexcept;
        void set_phase(const SampleIdType& sample, PhaseCall phase);
        
        virtual const GenomicRegion& mapped_region() const noexcept = 0;
        
        virtual const Allele& reference() const noexcept = 0;
        
        void replace(char old_base, char replacement_base);
        
        virtual void replace(const Allele& old, Allele replacement) = 0;
        virtual void replace_uncalled_genotype_alleles(const Allele& replacement, char ignore) = 0;
        
        virtual bool parsimonise(char dummy_base) { return false; };
        virtual bool parsimonise(const ReferenceGenome& reference) { return false; };
        
        virtual void decorate(VcfRecord::Builder& record) const = 0;
        
        void set_dummy_model_posterior(double p) noexcept;
        boost::optional<double> dummy_model_posterior() const noexcept;
        
    protected:
        std::unordered_map<SampleIdType, GenotypeCall> genotype_calls_;
        
        Phred<double> quality_;
        
        boost::optional<double> dummy_model_posterior_;
        
    private:
        virtual void replace_called_alleles(const char old_base, const char replacement_base) = 0;
    };
    
    template <typename T>
    Call::Call(T&& genotype_calls, Phred<double> quality)
    :
    genotype_calls_ {std::begin(genotype_calls),
    std::end(genotype_calls)},
    quality_ {quality},
    dummy_model_posterior_ {}
    {}
    
    template <typename R>
    Call::PhaseCall::PhaseCall(R&& region, Phred<double> score)
    : region_ {std::forward<R>(region)}, score_ {score}
    {}
    
    template <typename G>
    Call::GenotypeCall::GenotypeCall(G&& genotype, Phred<double> posterior)
    : genotype {std::forward<G>(genotype)}, posterior {posterior}, phase {}
    {}
    
    template <typename G, typename P>
    Call::GenotypeCall::GenotypeCall(G&& genotype, Phred<double> posterior, P&& phase)
    : genotype {std::forward<G>(genotype)}, posterior {posterior}, phase {std::forward<P>(phase)}
    {}
} // namespace Octopus

#endif /* call_hpp */
