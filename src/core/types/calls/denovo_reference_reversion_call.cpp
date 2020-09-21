// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "denovo_reference_reversion_call.hpp"

namespace octopus {

namespace {

class DummyGenerator
{
public:
    DummyGenerator() = delete;
    
    DummyGenerator(const char dummy) : dummy_ {dummy} {}
    
    char operator()(const GenomicRegion& region) const
    {
        return dummy_;
    }
private:
    char dummy_;
};

} // namespace

bool DenovoReferenceReversionCall::parsimonise(const char dummy_base)
{
    if (is_sequence_empty(variant_.ref_allele())) {
        const auto old_allele = variant_.ref_allele();
        variant_ = make_parsimonious(variant_, DummyGenerator {dummy_base});
        for (auto& p : genotype_calls_) {
            Genotype<Allele>& genotype {p.second.genotype};
            Genotype<Allele> parsimonised_genotype {genotype.ploidy()};
            for (const Allele& allele : genotype) {
                if (allele == old_allele) {
                    parsimonised_genotype.emplace(variant_.ref_allele());
                } else {
                    auto old_sequence = allele.sequence();
                    old_sequence.insert(std::begin(old_sequence), dummy_base);
                    Allele new_allele {octopus::mapped_region(variant_), std::move(old_sequence)};
                    parsimonised_genotype.emplace(std::move(new_allele));
                }
            }
            genotype = std::move(parsimonised_genotype);
        }
        return true;
    } else {
        return false;
    }
}

void DenovoReferenceReversionCall::decorate(VcfRecord::Builder& record) const
{
    DenovoCall::decorate(record);
    record.set_reference_reversion();
}

std::unique_ptr<Call> DenovoReferenceReversionCall::do_clone() const
{
    return std::make_unique<DenovoReferenceReversionCall>(*this);
}

} // namespace octopus
