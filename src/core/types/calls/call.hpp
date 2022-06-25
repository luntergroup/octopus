// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef call_hpp
#define call_hpp

#include <vector>
#include <unordered_map>
#include <utility>
#include <iterator>

#include <boost/optional.hpp>
#include <boost/container/flat_map.hpp>

#include "config/common.hpp"
#include "basics/phred.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "core/types/allele.hpp"
#include "core/types/genotype.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/variant/vcf_record.hpp"

namespace octopus {

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
        boost::optional<Phred<double>> model_posterior;
    };
    
    Call() = delete;
    
    explicit Call(Phred<double> quality);
    
    template <typename T> explicit Call(T&& genotype_calls, Phred<double> quality);
    
    Call(const Call&)            = default;
    Call& operator=(const Call&) = default;
    Call(Call&&)                 = default;
    Call& operator=(Call&&)      = default;
    
    virtual ~Call() = default;
    
    std::unique_ptr<Call> clone() const;
    
    Phred<double> quality() const noexcept;
    void set_quality(Phred<double> new_quality) noexcept;
    
    GenotypeCall& get_genotype_call(const SampleName& sample);
    const GenotypeCall& get_genotype_call(const SampleName& sample) const;
    
    bool is_phased(const SampleName& sample) const;
    bool all_phased() const noexcept;
    void set_phase(const SampleName& sample, PhaseCall phase);
    
    virtual const GenomicRegion& mapped_region() const noexcept = 0;
    
    virtual const Allele& reference() const noexcept = 0;
    
    virtual bool is_represented(const Allele& allele) const noexcept = 0;
    
    void replace(char old_base, char replacement_base);
    
    virtual void replace(const Allele& old, Allele replacement) = 0;
    virtual void replace_uncalled_genotype_alleles(const Allele& replacement, char ignore) = 0;
    
    virtual bool parsimonise(char dummy_base) { return false; };
    virtual bool parsimonise(const ReferenceGenome& reference) { return false; };
    
    virtual void decorate(VcfRecord::Builder& record) const = 0;
    
    virtual bool requires_model_evaluation() const noexcept { return false; };
    
    void set_model_posterior(Phred<double> p) noexcept;
    void set_model_posterior(const SampleName& sample, Phred<double> p) noexcept;
    
    boost::optional<Phred<double>> model_posterior() const noexcept;
    boost::optional<Phred<double>> model_posterior(const SampleName& sample) const noexcept;
    
    void reorder_genotype(const SampleName& sample, const std::vector<unsigned>& order);
    
protected:
    boost::container::flat_map<SampleName, GenotypeCall> genotype_calls_;
    Phred<double> quality_;
    boost::optional<Phred<double>> model_posterior_;
    
private:
    virtual std::unique_ptr<Call> do_clone() const = 0;
    virtual void replace_called_alleles(const char old_base, const char replacement_base) = 0;
    virtual void reorder_genotype_fields(const SampleName& sample, const std::vector<unsigned>& order) {};
};

template <typename T>
Call::Call(T&& genotype_calls, Phred<double> quality)
: genotype_calls_ {std::begin(genotype_calls), std::end(genotype_calls)}
, quality_ {quality}
, model_posterior_ {}
{}

template <typename R>
Call::PhaseCall::PhaseCall(R&& region, Phred<double> score)
: region_ {std::forward<R>(region)}
, score_ {score}
{}

template <typename G>
Call::GenotypeCall::GenotypeCall(G&& genotype, Phred<double> posterior)
: genotype {std::forward<G>(genotype)}
, posterior {posterior}
, phase {}
{}

template <typename G, typename P>
Call::GenotypeCall::GenotypeCall(G&& genotype, Phred<double> posterior, P&& phase)
: genotype {std::forward<G>(genotype)}
, posterior {posterior}
, phase {std::forward<P>(phase)}
{}

} // namespace octopus

#endif
