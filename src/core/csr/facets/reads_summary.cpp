// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "reads_summary.hpp"

#include <utility>

#include "utils/read_duplicates.hpp"

namespace octopus { namespace csr {

const std::string ReadsSummary::name_ {"ReadsSummary"};

auto copy_duplicates(const ReadContainer& reads)
{
    const auto duplicate_itrs = find_duplicate_reads(std::cbegin(reads), std::cend(reads));
    std::vector<Facet::ReadsSummary::DuplicateReadSet> result {};
    result.reserve(duplicate_itrs.size());
    for (const auto& itrs : duplicate_itrs) {
        assert(!itrs.empty());
        Facet::ReadsSummary::DuplicateReadSet dups {};
        dups.reads.reserve(itrs.size());
        std::transform(std::cbegin(itrs), std::cend(itrs), std::back_inserter(dups.reads), [] (auto itr) { return *itr; });
        result.push_back(std::move(dups));
    }
    std::sort(std::begin(result), std::end(result));
    return result;
}

ReadsSummary::ReadsSummary(const ReadMap& reads)
{
    result_.reserve(reads.size());
    for (const auto& p : reads) {
        result_[p.first].duplicates = copy_duplicates(p.second);
    }
}

Facet::ResultType ReadsSummary::do_get() const
{
    return std::cref(result_);
}

} // namespace csr
} // namespace octopus
