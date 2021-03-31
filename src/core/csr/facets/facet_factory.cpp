// Copyright (c) 2015-2020 Daniel Cooke
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
#include "repeat_context.hpp"
#include "samples.hpp"
#include "genotypes.hpp"
#include "alleles.hpp"
#include "ploidies.hpp"
#include "pedigree.hpp"
#include "reads_summary.hpp"

namespace octopus { namespace csr {

FacetFactory::FacetFactory(VcfHeader input_header)
: input_header_ {std::move(input_header)}
, samples_ {input_header_.samples()}
, reference_ {}
, read_pipe_ {}
, ploidies_ {}
, pedigree_ {}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory::FacetFactory(VcfHeader input_header,
                           const ReferenceGenome& reference,
                           BufferedReadPipe read_pipe,
                           PloidyMap ploidies,
                           HaplotypeLikelihoodModel likelihood_model)
: input_header_ {std::move(input_header)}
, samples_ {input_header_.samples()}
, reference_ {reference}
, read_pipe_ {std::move(read_pipe)}
, ploidies_ {std::move(ploidies)}
, pedigree_ {}
, likelihood_model_ {std::move(likelihood_model)}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory::FacetFactory(VcfHeader input_header,
                           const ReferenceGenome& reference,
                           BufferedReadPipe read_pipe,
                           PloidyMap ploidies,
                           HaplotypeLikelihoodModel likelihood_model,
                           octopus::Pedigree pedigree)
: input_header_ {std::move(input_header)}
, samples_ {input_header_.samples()}
, reference_ {reference}
, read_pipe_ {std::move(read_pipe)}
, ploidies_ {std::move(ploidies)}
, pedigree_ {std::move(pedigree)}
, likelihood_model_ {std::move(likelihood_model)}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory::FacetFactory(FacetFactory&& other)
: input_header_ {std::move(other.input_header_)}
, samples_ {std::move(other.samples_)}
, reference_ {std::move(other.reference_)}
, read_pipe_ {std::move(other.read_pipe_)}
, ploidies_ {std::move(other.ploidies_)}
, pedigree_ {std::move(other.pedigree_)}
, likelihood_model_ {std::move(other.likelihood_model_)}
, facet_makers_ {}
{
    setup_facet_makers();
}

FacetFactory& FacetFactory::operator=(FacetFactory&& other)
{
    using std::swap;
    swap(input_header_, other.input_header_);
    swap(samples_, other.samples_);
    swap(reference_, other.reference_);
    swap(read_pipe_, other.read_pipe_);
    swap(ploidies_, other.ploidies_);
    swap(pedigree_, other.pedigree_);
    swap(likelihood_model_, other.likelihood_model_);
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
    check_requirements(name);
    const auto block_data = make_block_data({name}, block);
    return make(name, block_data);
}

FacetFactory::FacetBlock FacetFactory::make(const std::vector<std::string>& names, const CallBlock& block) const
{
    if (names.empty()) return {};
    check_requirements(names);
    const auto block_data = make_block_data(names, block);
    return make(names, block_data);
}

namespace {

template <typename Facet>
decltype(auto) name() noexcept
{
    return Facet().name();
}

bool requires_reference(const std::string& facet) noexcept
{
    const static std::array<std::string, 3> read_facets{name<ReferenceContext>(), name<RepeatContext>(), name<ReadAssignments>()};
    return std::find(std::cbegin(read_facets), std::cend(read_facets), facet) != std::cend(read_facets);
}

bool requires_reference(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [](const auto& facet) { return requires_reference(facet); });
}

bool requires_reads(const std::string& facet) noexcept
{
    const static std::array<std::string, 3> read_facets{name<OverlappingReads>(), name<ReadsSummary>(), name<ReadAssignments>()};
    return std::find(std::cbegin(read_facets), std::cend(read_facets), facet) != std::cend(read_facets);
}

bool requires_reads(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [](const auto& facet) { return requires_reads(facet); });
}

bool requires_genotypes(const std::string& facet) noexcept
{
    const static std::array<std::string, 2> genotype_facets {name<Genotypes>(), name<ReadAssignments>()};
    return std::find(std::cbegin(genotype_facets), std::cend(genotype_facets), facet) != std::cend(genotype_facets);
}

bool requires_genotypes(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [](const auto& facet) { return requires_genotypes(facet); });
}

bool requires_ploidies(const std::string& facet) noexcept
{
    const static std::array<std::string, 1> genotype_facets{name<Ploidies>()};
    return std::find(std::cbegin(genotype_facets), std::cend(genotype_facets), facet) != std::cend(genotype_facets);
}

bool requires_ploidies(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [](const auto& facet) { return requires_ploidies(facet); });
}

bool requires_pedigree(const std::string& facet) noexcept
{
    const static std::array<std::string, 1> genotype_facets{name<Pedigree>()};
    return std::find(std::cbegin(genotype_facets), std::cend(genotype_facets), facet) != std::cend(genotype_facets);
}

bool requires_pedigree(const std::vector<std::string>& facets) noexcept
{
    return std::any_of(std::cbegin(facets), std::cend(facets), [](const auto& facet) { return requires_pedigree(facet); });
}

} // namespace

class BadFacetFactoryRequest : public ProgramError
{
    std::string facet_;
    std::string do_where() const override { return "FacetFactory"; }
    std::string do_why() const override
    {
        return "Could not make facet " + facet_ + " due to an unmet requirement";
    }
    std::string do_help() const override
    {
        return "submit an error report";
    }
public:
    BadFacetFactoryRequest(std::string facet) : facet_ {std::move(facet)} {}
};

std::vector<FacetFactory::FacetBlock>
FacetFactory::make(const std::vector<std::string>& names,
                   const std::vector<CallBlock>& blocks,
                   ThreadPool& workers) const
{
    if (blocks.empty()) return {};
    check_requirements(names);
    std::vector<FacetBlock> result {};
    result.reserve(blocks.size());
    if (blocks.size() > 1 && !workers.empty()) {
        std::vector<std::future<FacetBlock>> futures {};
        futures.reserve(blocks.size());
        const auto fetch_reads = requires_reads(names);
        const auto fetch_genotypes = requires_genotypes(names);
        for (const auto& block : blocks) {
            // It's faster to fetch reads serially from left to right, so do this outside the thread pool
            BlockData data {};
            data.calls = std::addressof(block);
            if (!block.empty()) {
                data.region = encompassing_region(block);
                if (fetch_reads) {
                    data.reads = read_pipe_->fetch_reads(*data.region);
                }
            }
            futures.push_back(workers.push([this, &names, data {std::move(data)}, &block, fetch_genotypes] () mutable {
                if (fetch_genotypes) {
                    data.genotypes = extract_genotypes(block, samples_, *reference_);
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
    facet_makers_[name<ReadsSummary>()] = [] (const BlockData& block) -> FacetWrapper
    {
        assert(block.reads);
        return {std::make_unique<ReadsSummary>(*block.reads)};
    };
    facet_makers_[name<ReadAssignments>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        assert(block.reads && block.genotypes);
        if (likelihood_model_) {
            return {std::make_unique<ReadAssignments>(*reference_, *block.genotypes, *block.reads, *block.calls, *likelihood_model_)};
        } else {
            return {std::make_unique<ReadAssignments>(*reference_, *block.genotypes, *block.reads, *block.calls)};
        }
    };
    facet_makers_[name<ReferenceContext>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        if (block.region) {
            constexpr GenomicRegion::Size context_size {50};
            return {std::make_unique<ReferenceContext>(*reference_, expand(*block.region, context_size))};
        } else {
            return {nullptr};
        }
    };
    facet_makers_[name<RepeatContext>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        if (block.region) {
            constexpr GenomicRegion::Size context_size {50};
            return {std::make_unique<RepeatContext>(*reference_, expand(*block.region, context_size))};
        } else {
            return {nullptr};
        }
    };
    facet_makers_[name<Samples>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        return {std::make_unique<Samples>(this->samples_)};
    };
    facet_makers_[name<Genotypes>()] = [] (const BlockData& block) -> FacetWrapper
    {
        assert(block.genotypes);
        return {std::make_unique<Genotypes>(*block.genotypes)};
    };
    facet_makers_[name<Alleles>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        assert(block.calls);
        return {std::make_unique<Alleles>(this->samples_, *block.calls)};
    };
    facet_makers_[name<Ploidies>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        return {std::make_unique<Ploidies>(*ploidies_, *block.region, input_header_.samples())};
    };
    facet_makers_[name<Pedigree>()] = [this] (const BlockData& block) -> FacetWrapper
    {
        return {std::make_unique<Pedigree>(*pedigree_)};
    };
}

void FacetFactory::check_requirements(const std::string& name) const
{
    if (!read_pipe_ && requires_reads(name)) {
        throw BadFacetFactoryRequest {name};
    }
    if (!reference_ && requires_reference(name)) {
        throw BadFacetFactoryRequest {name};
    }
    if (!ploidies_ && requires_ploidies(name)) {
        throw BadFacetFactoryRequest {name};
    }
    if (!pedigree_ && requires_pedigree(name)) {
        throw BadFacetFactoryRequest {name};
    }
}

void FacetFactory::check_requirements(const std::vector<std::string>& names) const
{
    for (const auto& name : names) {
        check_requirements(name);
    }
}

FacetWrapper FacetFactory::make(const std::string& name, const BlockData& block) const
{
    if (facet_makers_.count(name) == 1) {
        return facet_makers_.at(name)(block);
    } else {
        throw UnknownFacet {name};
    }
}

FacetFactory::FacetBlock 
FacetFactory::make(const std::vector<std::string>& names, const BlockData& block) const
{
    FacetBlock result {};
    result.reserve(names.size());
    for (const auto& name : names) {
        result.push_back(make(name, block));
    }
    return result;
}

FacetFactory::FacetBlock
FacetFactory::make(const std::vector<std::string>& names, const BlockData& block,
                   ThreadPool& workers) const
{
    if (workers.empty() || names.size() < 2) return make(names, block);
    std::vector<std::future<FacetWrapper>> futures {};
    futures.reserve(names.size());
    for (const auto& name : names) {
        futures.push_back(workers.push([&] () mutable { return this->make(name, block); }));
    }
    FacetBlock result {};
    result.reserve(names.size());
    for (auto& fut : futures) {
        result.push_back(fut.get());
    }
    return result;
}

FacetFactory::BlockData FacetFactory::make_block_data(const std::vector<std::string>& names, const CallBlock& block) const
{
    BlockData result {};
    result.calls = std::addressof(block);
    assert(!names.empty());
    if (!block.empty()) {
        result.region = encompassing_region(block);
        if (requires_reads(names)) {
            result.reads = read_pipe_->fetch_reads(*result.region);
        }
        if (requires_genotypes(names)) {
            result.genotypes = extract_genotypes(block, samples_, *reference_);
        }
    }
    return result;
}

} // namespace csr
} // namespace octopus
