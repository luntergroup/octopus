// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cnv_call.hpp"

#include "utils/string_utils.hpp"
#include "utils/reorder.hpp"

namespace octopus {

void CNVCall::decorate(VcfRecord::Builder& record) const
{
    if (posterior_) {
        record.set_info("PP", utils::to_string(posterior_->score()));
    }
    if (!somatic_haplotypes_.empty()) {
        record.add_format({"HSS"});
        for (const auto& p : somatic_haplotypes_) {
            const auto ploidy = this->genotype_calls_.at(p.first).genotype.ploidy();
            assert(ploidy == p.second.size());
            std::vector<std::string> somatic_status(ploidy);
            std::transform(std::cbegin(p.second), std::cend(p.second), std::begin(somatic_status),
                           [] (bool status) { return std::to_string(status); });
            record.set_format(p.first, "HSS", std::move(somatic_status));
        }
    }
}

std::unique_ptr<Call> CNVCall::do_clone() const
{
    return std::make_unique<CNVCall>(*this);
}

void CNVCall::reorder_genotype_fields(const SampleName& sample, const std::vector<unsigned>& order)
{
    auto& somatics = somatic_haplotypes_.at(sample);
    assert(somatics.size() == order.size());
    reorder(std::cbegin(order), std::cend(order), std::begin(somatics));
}

} // namespace octopus
