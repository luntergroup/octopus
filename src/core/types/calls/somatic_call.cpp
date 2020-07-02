// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_call.hpp"

#include "utils/concat.hpp"
#include "utils/string_utils.hpp"
#include "utils/reorder.hpp"

namespace octopus {

SomaticCall::SomaticCall(Variant variant,
                         const CancerGenotype<Allele>& genotype_call,
                         Phred<double> quality,
                         Phred<double> genotype_posterior,
                         GenotypeStatsMap stats,
                         boost::optional<Phred<double>> classification_posterior)
: VariantCall {std::move(variant), decltype(genotype_calls_) {}, quality, classification_posterior}
, genotype_stats_ {}
{
    if (variant_.ref_allele() == variant_.alt_allele()) {
        Allele::NucleotideSequence missing_sequence(ref_sequence_size(variant_), 'N');
        using octopus::mapped_region;
        variant_ = Variant {Allele {mapped_region(variant_), std::move(missing_sequence)},variant_.alt_allele()};
    }
    genotype_calls_.reserve(stats.size()); // num samples
    genotype_stats_.reserve(stats.size());
    for (auto& p : stats) {
        auto& sample_stats = p.second;
        genotype_calls_.emplace(p.first, GenotypeCall {sample_stats.somatic.empty() ? genotype_call.germline() : demote(genotype_call), genotype_posterior});
        genotype_stats_.emplace(p.first, concat(std::move(sample_stats.germline), std::move(sample_stats.somatic)));
    }
}

namespace {

static std::string to_string_sf(const double val, const int sf = 2)
{
    return utils::strip_leading_zeroes(utils::to_string(val, sf, utils::PrecisionRule::sf));
}

} // namespace

void SomaticCall::decorate(VcfRecord::Builder& record) const
{
    record.set_somatic();
    if (posterior_) {
        record.set_info("PP", utils::to_string(posterior_->score()));
    }
    record.add_format({"HPC", "MAP_HF", "HF_CR"});
    for (const auto& p : genotype_stats_) {
        std::vector<std::string> hpcs {}, map_hfs {}, hf_crs {};
        const auto ploidy = this->genotype_calls_.at(p.first).genotype.ploidy();
        assert(ploidy == p.second.size());
        hpcs.reserve(ploidy);
        map_hfs.reserve(ploidy);
        hf_crs.reserve(2 * ploidy);
        for (const auto& stats : p.second) {
            hpcs.push_back(to_string_sf(stats.pseudo_count));
            map_hfs.push_back(to_string_sf(stats.map_vaf));
            hf_crs.insert(std::cend(hf_crs), {to_string_sf(stats.vaf_credible_region.first), to_string_sf(stats.vaf_credible_region.second)});
        }
        record.set_format(p.first, "HPC", std::move(hpcs));
        record.set_format(p.first, "MAP_HF", std::move(map_hfs));
        record.set_format(p.first, "HF_CR", std::move(hf_crs));
    }
}

std::unique_ptr<Call> SomaticCall::do_clone() const
{
    return std::make_unique<SomaticCall>(*this);
}

void SomaticCall::reorder_genotype_fields(const SampleName& sample, const std::vector<unsigned>& order)
{
    auto& stats = genotype_stats_.at(sample);
    assert(stats.size() == order.size());
    reorder(std::cbegin(order), std::cend(order), std::begin(stats));
}

} // namespace octopus
