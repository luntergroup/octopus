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
    alleles_.reserve(calls.size());
    for (const auto& call : calls) {
        auto& call_map = alleles_[mapped_region(call)];
        for (const auto& sample : samples) {
            call_map[sample] = get_called_alleles(call, sample);
        }
    }
}

Facet::ResultType Alleles::do_get() const
{
    return std::cref(alleles_);
}

std::vector<Allele> get_all(const Facet::AlleleMap& alleles, const VcfRecord& call, const SampleName& sample)
{
    return alleles.at(mapped_region(call)).at(sample).first;
}

std::vector<Allele> get_alt(const Facet::AlleleMap& alleles, const VcfRecord& call, const SampleName& sample)
{
    const auto& p = alleles.at(mapped_region(call)).at(sample);
    auto result = p.first;
    if (p.second) result.erase(std::cbegin(result));
    return result;
}

std::vector<Allele> get_all_unique(const Facet::AlleleMap& alleles, const VcfRecord& call, const std::vector<SampleName>& samples)
{
    std::vector<Allele> result {};
    for (const auto& sample : samples) {
        utils::append(get_all(alleles, call, sample), result);
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

} // namespace csr
} // namespace octopus
