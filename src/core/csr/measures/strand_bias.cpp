// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "strand_bias.hpp"

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "basics/aligned_read.hpp"
#include "utils/maths.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

std::unique_ptr<Measure> StrandBias::do_clone() const
{
    return std::make_unique<StrandBias>(*this);
}

bool is_forward(const AlignedRead& read) noexcept
{
    return read.direction() == AlignedRead::Direction::forward;
}

template <typename Container>
auto count_directions(const Container& reads)
{
    auto n_fwd = static_cast<std::size_t>(std::count_if(std::cbegin(reads), std::cend(reads),
                                                        [] (const auto& read) { return is_forward(read); }));
    return std::make_pair(n_fwd, reads.size() - n_fwd);
}

Measure::ResultType StrandBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto assignments = boost::get<ReadAssignments::ResultType>(facets.at("ReadAssignments").get());
    for (const auto& p : assignments) {
        if (call.is_heterozygous(p.first)) {
            std::vector<std::pair<unsigned, unsigned>> direction_counts {};
            direction_counts.reserve(p.second.size());
            for (const auto& h : p.second) {
                direction_counts.push_back(count_directions(h.second));
            }
        }
    }
    return boost::none;
}

std::string StrandBias::do_name() const
{
    return "SB";
}

std::vector<std::string> StrandBias::do_requirements() const
{
    return {"ReadAssignments"};
}

} // namespace csr
} // namespace octopus