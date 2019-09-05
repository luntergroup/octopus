// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "alleles.hpp"

#include <iterator>
#include <algorithm>

#include "utils/genotype_reader.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/append.hpp"

namespace octopus { namespace csr {

const std::string Alleles::name_ {"Alleles"};

template <typename T>
static void insert(std::vector<T>&& src, MappableFlatSet<T>& dst)
{
    dst.insert(std::make_move_iterator(std::begin(src)), std::make_move_iterator(std::end(src)));
}

Alleles::Alleles(const std::vector<VcfRecord::SampleName>& samples, const std::vector<VcfRecord>& calls)
{
    alleles_.reserve(samples.size());
    for (const auto& sample : samples) {
        for (const auto& call : calls) {
            insert(get_called_alleles(call, sample).first, alleles_[sample]);
        }
    }
}

Facet::ResultType Alleles::do_get() const
{
    return std::cref(alleles_);
}

std::vector<Allele> copy_overlapped(const MappableFlatSet<Allele>& alleles, const VcfRecord& call)
{
    const auto overlapped = overlap_range(alleles, call);
    return {std::cbegin(overlapped), std::cend(overlapped)};
}

std::vector<Allele> copy_unique_overlapped(const Facet::AlleleMap& alleles, const VcfRecord& call, const std::vector<SampleName>& samples)
{
    std::vector<Allele> result {};
    result.reserve(2 * alleles.size());
    for (const auto& sample : samples) {
        utils::append(copy_overlapped(alleles.at(sample), call), result);
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

} // namespace csr
} // namespace octopus
