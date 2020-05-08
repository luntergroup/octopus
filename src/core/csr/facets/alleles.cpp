// Copyright (c) 2015-2020 Daniel Cooke
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
    alleles_.reserve(calls.size());
    for (const auto& call : calls) {
        auto& call_map = alleles_[mapped_region(call)];
        for (const auto& sample : samples) {
            call_map[sample] = get_resolved_alleles(call, sample);
        }
    }
}

Facet::ResultType Alleles::do_get() const
{
    return std::cref(alleles_);
}

std::vector<Allele> get_called(const Facet::AlleleMap& alleles, const VcfRecord& call, const SampleName& sample)
{
    std::vector<Allele> result {};
    result.reserve(call.alt().size() + 1);
    for (const auto& allele : get(alleles, call, sample)) {
        if (allele) result.push_back(*allele);
    }
    return result;
}

std::vector<Allele> get_unique_called(const Facet::AlleleMap& alleles, const VcfRecord& call, const std::vector<SampleName>& samples)
{
    std::vector<Allele> result {};
    result.reserve(call.alt().size() + 1);
    for (const auto& sample : samples) {
        for (const auto& allele : get(alleles, call, sample)) {
            if (allele) result.push_back(*allele);
        }
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

} // namespace csr
} // namespace octopus
