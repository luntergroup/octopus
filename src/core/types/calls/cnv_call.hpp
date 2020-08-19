// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef cnv_call_hpp
#define cnv_call_hpp

#include <utility>
#include <unordered_map>
#include <vector>

#include <boost/optional.hpp>
#include <boost/container/flat_map.hpp>

#include "config/common.hpp"
#include "core/types/allele.hpp"
#include "core/types/variant.hpp"
#include "germline_variant_call.hpp"

namespace octopus {

class CNVCall : public GermlineVariantCall
{
public:    
    using SomaticHaplotypeVector = std::vector<bool>;
    using SomaticHaplotypeMap = std::unordered_map<SampleName, SomaticHaplotypeVector>;

    CNVCall() = delete;
    
    template <typename V, typename T>
    CNVCall(V&& variant, T&& genotype_calls, Phred<double> quality,
            boost::optional<Phred<double>> posterior = boost::none,
            boost::optional<SomaticHaplotypeMap> somatic_haplotypes = boost::none);
    
    CNVCall(const CNVCall&)            = default;
    CNVCall& operator=(const CNVCall&) = default;
    CNVCall(CNVCall&&)                 = default;
    CNVCall& operator=(CNVCall&&)      = default;
    
    virtual ~CNVCall() = default;
    virtual void decorate(VcfRecord::Builder& record) const override;
    
protected:
    boost::container::flat_map<SampleName, SomaticHaplotypeVector> somatic_haplotypes_;
    
private:
    virtual std::unique_ptr<Call> do_clone() const override;
    virtual void reorder_genotype_fields(const SampleName& sample, const std::vector<unsigned>& order) override;
};

template <typename V, typename T>
CNVCall::CNVCall(V&& variant, T&& genotype_calls, Phred<double> quality,
                 boost::optional<Phred<double>> posterior,
                 boost::optional<SomaticHaplotypeMap> somatic_haplotypes)
: GermlineVariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality, posterior}
, somatic_haplotypes_ {}
{
    if (somatic_haplotypes) {
        somatic_haplotypes_.reserve(somatic_haplotypes->size());
        for (auto& p : *somatic_haplotypes) {
            somatic_haplotypes_.emplace(p.first, std::move(p.second));
        }
    }
}

} // namespace octopus

#endif
