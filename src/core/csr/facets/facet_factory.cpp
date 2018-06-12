// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "facet_factory.hpp"

#include <utility>
#include <memory>
#include <algorithm>
#include <iterator>
#include <array>
#include <future>
#include <cassert>

#include "exceptions/program_error.hpp"
#include "overlapping_reads.hpp"
#include "read_assignments.hpp"
#include "reference_context.hpp"
#include "samples.hpp"
#include "genotypes.hpp"
#include "ploidies.hpp"

namespace octopus { namespace csr {

FacetFactory::FacetFactory(const ReferenceGenome& reference, BufferedReadPipe read_pipe, VcfHeader input_header, PloidyMap ploidies)
: reference_ {reference}
, read_pipe_ {std::move(read_pipe)}
, input_header_ {std::move(input_header)}
, ploidies_ {std::move(ploidies)}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory::FacetFactory(FacetFactory&& other)
: reference_ {std::move(other.reference_)}
, read_pipe_ {std::move(other.read_pipe_)}
, input_header_ {std::move(other.input_header_)}
, ploidies_ {std::move(other.ploidies_)}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory& FacetFactory::operator=(FacetFactory&& other)
{
    using std::swap;
    swap(reference_, other.reference_);
    swap(read_pipe_, other.read_pipe_);
    swap(ploidies_, other.ploidies_);
    swap(input_header_, other.input_header_);
    setup_facet_makers();
    return *this;
}

class UnknownFacet : public ProgramError
{
    std::string do_where() const override { return "FacetFactory::make"; }
    std::string do_why() const override
    {
        return std::string {"Trying to make unknown Facet "} + name_;
    }

public:
    UnknownFacet(std::string name) : name_ {std::move(name)} {}
    std::string name_;
};

FacetWrapper FacetFactory::make(const std::string& name, const CallBlock& block) const
{
    const auto block_data = make_block_data({name}, block);
    return make(name, block_data);
}

FacetFactory::FacetBlock FacetFactory::make(const std::vector<std::string>& names, const CallBlock& block) const
{
    if (names.empty()) return {};
    const auto block_data = make_block_data(names, block);
    return make(names, block_data);
}

namespace {

template <typename Facet>
decltype(auto) name() noexcept
{
    return Facet().name();
}

bool requires_reads(const std::string& facet) noexcept
{
    const static std::array<std::string, 2> read_facets{name<OverlappingReads>(), name<ReadAssignments>()};
    return std::find(std::cbegin(read_facets), std::cend(read_facets), facet) != std::cend(read_facets);
}

bool requires_reads(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [](const auto& facet) { return requires_reads(facet); });
}

bool requires_genotypes(const std::string& facet) noexcept
{
    const static std::array<std::string, 1> genotype_facets{name<ReadAssignments>()};
    return std::find(std::cbegin(genotype_facets), std::cend(genotype_facets), facet) != std::cend(genotype_facets);
}

bool requires_genotypes(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets),
                       [](const auto& facet) { return requires_genotypes(facet); });
}

} // namespace

std::vector<FacetFactory::FacetBlock> FacetFactory::make(const std::vector<std::string>& names,
                                                         const std::vector<CallBlock>& blocks,
                                                         ThreadPool& workers) const
{
    if (blocks.empty()) return {};
    std::vector<FacetBlock> result {};
    result.reserve(blocks.size());
    if (blocks.size() > 1 && !workers.empty()) {
        std::vector<std::future<FacetBlock>> futures {};
        futures.reserve(blocks.size());
        const auto fetch_reads = requires_reads(names);
        const auto fetch_genotypes = requires_genotypes(names);
        const auto& samples = read_pipe_.source().samples();
        for (const auto& block : blocks) {
            // It's faster to fetch reads serially from left to right, so do this outside the thread pool
            BlockData data {};
            if (!block.empty()) {
                data.region = encompassing_region(block);
                if (fetch_reads) {
                    data.reads = read_pipe_.fetch_reads(*data.region);
                }
            }
            futures.push_back(workers.push([this, &names, data {std::move(data)}, &block, &samples, fetch_genotypes] () mutable {
                if (fetch_genotypes) {
                    data.genotypes = extract_genotypes(block, samples, reference_);
                }
                return this->make(names, data);
            }));
        }
        for (auto& fut : futures) {
            result.push_back(fut.get());
        }
    } else {
        for (const auto& block : blocks) {
            const auto data = make_block_data(names, block);
            result.push_back(make(names, data));
        }
    }
    return result;
}

// private methods

void FacetFactory::setup_facet_makers()
{
    facet_makers_[name<OverlappingReads>()] = [] (const BlockData& block) -> FacetWrapper
    {
        assert(block.reads);
        return {std::make_unique<OverlappingReads>(*block.reads)};
    };
    facet_makers_[name<ReadAssignments>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        assert(block.reads && block.genotypes);
        return {std::make_unique<ReadAssignments>(reference_, *block.genotypes, *block.reads)};
    };
    facet_makers_[name<ReferenceContext>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        if (block.region) {
            constexpr GenomicRegion::Size context_size {50};
            return {std::make_unique<ReferenceContext>(reference_, expand(*block.region, context_size))};
        } else {
            return {nullptr};
        }
    };
    facet_makers_[name<Samples>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        return {std::make_unique<Samples>(input_header_.samples())};
    };
    facet_makers_[name<Genotypes>()] = [] (const BlockData& block) -> FacetWrapper
    {
        assert(block.genotypes);
        return {std::make_unique<Genotypes>(*block.genotypes)};
    };
    facet_makers_[name<Ploidies>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        return {std::make_unique<Ploidies>(ploidies_, *block.region, input_header_.samples())};
    };
}

FacetWrapper FacetFactory::make(const std::string& name, const BlockData& block) const
{
    if (facet_makers_.count(name) == 1) {
        return facet_makers_.at(name)(block);
    } else {
        throw UnknownFacet {name};
    }
}

FacetFactory::FacetBlock FacetFactory::make(const std::vector<std::string>& names, const BlockData& block) const
{
    FacetBlock result {};
    result.reserve(names.size());
    for (const auto& name : names) {
        result.push_back(make(name, block));
    }
    return result;
}

FacetFactory::BlockData FacetFactory::make_block_data(const std::vector<std::string>& names, const CallBlock& block) const
{
    BlockData result {};
    assert(!names.empty());
    if (!block.empty()) {
        result.region = encompassing_region(block);
        if (requires_reads(names)) {
            result.reads = read_pipe_.fetch_reads(*result.region);
        }
        if (requires_genotypes(names)) {
            result.genotypes = extract_genotypes(block, read_pipe_.source().samples(), reference_);
        }
    }
    return result;
}

} // namespace csr
} // namespace octopus
