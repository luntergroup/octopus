// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cell_variant_call.hpp"

#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>
#include <iostream>

#include "utils/string_utils.hpp"

namespace octopus {

void CellVariantCall::decorate(VcfRecord::Builder& record) const
{
    if (is_somatic()) record.set_somatic();
    record.set_info("PPP", utils::to_string(phylogeny_summary_.map_posterior.score()));
    std::vector<std::string> size_posteriors(phylogeny_summary_.size_posteriors.size() - 1);
    std::transform(std::next(std::cbegin(phylogeny_summary_.size_posteriors)), std::cend(phylogeny_summary_.size_posteriors),
                   std::begin(size_posteriors), [] (auto p) { return std::to_string(static_cast<unsigned>(p.score())); });
    record.set_info("PSPP", std::move(size_posteriors));
    if (phylogeny_summary_.map.size() > 1) {
        std::ostringstream ss {};
        phylogeny_summary_.map.serialise(ss, [] (std::ostream& os, const auto& group) { os << group.id; });
        record.set_info("PY", ss.str());
        record.add_format("PNAP");
        for (const auto& p : phylogeny_summary_.sample_node_posteriors) {
            assert(p.second.size() == phylogeny_summary_.map.size());
            std::vector<std::string> posteriors {};
            posteriors.reserve(p.second.size());
            for (auto posterior : p.second) posteriors.push_back(utils::to_string(std::min(posterior.score(), 20.0)));
            record.set_format(p.first, "PNAP", std::move(posteriors));
        }
    }
}

std::unique_ptr<Call> CellVariantCall::do_clone() const
{
    return std::make_unique<CellVariantCall>(*this);
}

bool CellVariantCall::is_somatic() const
{
    const auto genotype_call_not_equal = [] (const auto& lhs, const auto& rhs) { return !are_equivalent(lhs.second.genotype, rhs.second.genotype); };
    return std::adjacent_find(std::cbegin(genotype_calls_), std::cend(genotype_calls_), genotype_call_not_equal) != std::cend(genotype_calls_);
}

} // namespace octopus
