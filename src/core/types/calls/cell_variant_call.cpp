// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cell_variant_call.hpp"

#include <iterator>
#include <algorithm>

#include "utils/string_utils.hpp"

namespace octopus {

void CellVariantCall::decorate(VcfRecord::Builder& record) const
{
    if (is_somatic()) {
        record.set_somatic();
    }
    record.set_info("PPP", utils::to_string(phylogeny_summary_.map_posterior.score()));
    std::vector<std::string> size_posteriors {};
    size_posteriors.reserve(phylogeny_summary_.size_posteriors.size());
    for (const auto p : phylogeny_summary_.size_posteriors) {
        size_posteriors.push_back(std::to_string(static_cast<unsigned>(p.score())));
    }
    record.set_info("PSPP", std::move(size_posteriors));
}

std::unique_ptr<Call> CellVariantCall::do_clone() const
{
    return std::make_unique<CellVariantCall>(*this);
}

bool CellVariantCall::is_somatic() const
{
    const auto genotype_call_not_equal = [] (const auto& lhs, const auto& rhs) { return lhs.second.genotype != rhs.second.genotype; };
    return std::adjacent_find(std::cbegin(genotype_calls_), std::cend(genotype_calls_), genotype_call_not_equal) != std::cend(genotype_calls_);
}

} // namespace octopus
