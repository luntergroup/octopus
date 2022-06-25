// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "randomiser.hpp"

#include <random>
#include <chrono>

#include "io/reference/reference_genome.hpp"
#include "basics/genomic_region.hpp"
#include "core/types/allele.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/sequence_utils.hpp"

namespace octopus { namespace coretools {

Randomiser::Randomiser(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
{}

std::unique_ptr<VariantGenerator> Randomiser::do_clone() const
{
    return std::make_unique<Randomiser>(*this);
}

void Randomiser::do_add_reads(const SampleName& sample, ReadVectorIterator first, ReadVectorIterator last)
{
    max_read_size_ = region_size(*largest_mappable(first, last));
}

void Randomiser::do_add_reads(const SampleName& sample, ReadFlatSetIterator first, ReadFlatSetIterator last)
{
    max_read_size_ = region_size(*largest_mappable(first, last));
}

std::vector<Variant> Randomiser::do_generate(const RegionSet& regions, OptionalThreadPool workers) const
{
    std::vector<Variant> result {};
    for (const auto& region : regions) {
        auto num_positions = region_size(region);
        std::vector<Variant> result {};
        if (num_positions == 0) return result;
        static const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
        static std::default_random_engine generator {static_cast<unsigned>(seed)};
        using T = Variant::MappingDomain::Size;
        std::uniform_int_distribution<T> uniform {0, std::min(num_positions, max_read_size_)};
        auto positions = decompose(region);
        for (auto p = uniform(generator); p < num_positions; p += max_read_size_) {
            auto position = positions[p];
            auto reference_allele = make_reference_allele(position, reference_);
            Allele mutation {position, utils::reverse_complement_copy(reference_allele.sequence())};
            result.emplace_back(reference_allele, mutation);
        }
    }
    return result;
}

std::string Randomiser::name() const
{
    return "Random";
}

} // namespace coretools
} // namespace octopus
