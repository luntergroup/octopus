// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "facet_factory.hpp"

#include <utility>
#include <memory>
#include <algorithm>
#include <cassert>

#include "exceptions/program_error.hpp"
#include "utils/parallel_transform.hpp"
#include "overlapping_reads.hpp"
#include "read_assignments.hpp"
#include "reference_context.hpp"
#include "samples.hpp"

namespace octopus { namespace csr {

FacetFactory::FacetFactory(const ReferenceGenome& reference, BufferedReadPipe read_pipe)
: reference_ {reference}
, read_pipe_ {std::move(read_pipe)}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory::FacetFactory(FacetFactory&& other)
: reference_ {std::move(other.reference_)}
, read_pipe_ {std::move(other.read_pipe_)}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory& FacetFactory::operator=(FacetFactory&& other)
{
    using std::swap;
    swap(reference_, other.reference_);
    swap(read_pipe_, other.read_pipe_);
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

std::vector<FacetFactory::FacetBlock> FacetFactory::make(const std::vector<std::string>& names,
                                                         const std::vector<CallBlock>& blocks,
                                                         ThreadPool& workers) const
{
    if (blocks.empty()) return {};
    std::vector<BlockData> data {};
    data.reserve(blocks.size());
    for (const auto& block : blocks) {
        data.push_back(make_block_data(names, block));
    }
    std::vector<FacetBlock> result {};
    result.reserve(blocks.size());
    if (blocks.size() > 1 && !workers.empty()) {
        transform(std::cbegin(data), std::cend(data), std::back_inserter(result),
                  [&] (const auto& block_data) { return make(names, block_data); },
                  workers);
    } else {
        std::transform(std::cbegin(data), std::cend(data), std::back_inserter(result),
                       [&] (const auto& block_data) { return make(names, block_data); });
    }
    return result;
}

// private methods

void FacetFactory::setup_facet_makers()
{
    facet_makers_["OverlappingReads"] = [] (const BlockData& block) -> FacetWrapper
    {
        assert(block.reads);
        return {std::make_unique<OverlappingReads>(*block.reads)};
    };
    facet_makers_["ReadAssignments"] = [this] (const BlockData& block) -> FacetWrapper
    {
        assert(block.reads && block.genotypes);
        return {std::make_unique<ReadAssignments>(reference_, *block.genotypes, *block.reads)};
    };
    facet_makers_["ReferenceContext"] = [this] (const BlockData& block) -> FacetWrapper
    {
        if (block.region) {
            constexpr GenomicRegion::Size context_size {50};
            return {std::make_unique<ReferenceContext>(reference_, expand(*block.region, context_size))};
        } else {
            return {nullptr};
        }
    };
    facet_makers_["Samples"] = [this] (const BlockData& block) -> FacetWrapper
    {
        return {std::make_unique<Samples>(read_pipe_.source().samples())};
    };
}

bool requires_reads(const std::string& facet) noexcept
{
    return facet == "OverlappingReads" || facet == "ReadAssignments";
}

bool requires_reads(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [] (const auto& facet) { return requires_reads(facet); });
}

bool requires_genotypes(const std::string& facet) noexcept
{
    return facet == "ReadAssignments";
}

bool requires_genotypes(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [] (const auto& facet) { return requires_genotypes(facet); });
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
