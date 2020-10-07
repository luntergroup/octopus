// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "polyclone_variant_call.hpp"

#include "utils/string_utils.hpp"
#include "utils/reorder.hpp"

namespace octopus {

namespace {

static std::string to_string_sf(const double val, const int sf = 2)
{
    return utils::strip_leading_zeroes(utils::to_string(val, sf, utils::PrecisionRule::sf));
}

} // namespace

void PolycloneVariantCall::decorate(VcfRecord::Builder& record) const
{
    if (posterior_) {
        record.set_info("PP", utils::to_string(posterior_->score()));
    }
    if (haplotype_frequency_stats_ && haplotype_frequency_stats_->size() > 1) {
        record.add_format({"HPC", "MAP_HF"});
        std::vector<std::string> hpcs {}, map_hfs {};
        const auto ploidy = std::cbegin(this->genotype_calls_)->second.genotype.ploidy();
        assert(ploidy == haplotype_frequency_stats_->size());
        hpcs.reserve(ploidy);
        map_hfs.reserve(ploidy);
        for (const auto& stats : *haplotype_frequency_stats_) {
            hpcs.push_back(to_string_sf(stats.pseudo_count));
            map_hfs.push_back(to_string_sf(stats.map));
        }
        const auto& sample = std::cbegin(this->genotype_calls_)->first;
        record.set_format(sample, "HPC", std::move(hpcs));
        record.set_format(sample, "MAP_HF", std::move(map_hfs));
    }
}

std::unique_ptr<Call> PolycloneVariantCall::do_clone() const
{
    return std::make_unique<PolycloneVariantCall>(*this);
}

void PolycloneVariantCall::reorder_genotype_fields(const SampleName& sample, const std::vector<unsigned>& order)
{
    if (haplotype_frequency_stats_) {
        assert(haplotype_frequency_stats_->size() == order.size());
        reorder(std::cbegin(order), std::cend(order), std::begin(*haplotype_frequency_stats_));
    }    
}

} // namespace octopus
