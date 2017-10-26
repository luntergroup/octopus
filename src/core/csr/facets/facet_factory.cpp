// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "facet_factory.hpp"

#include <utility>
#include <memory>
#include <algorithm>

#include "read_assignments.hpp"
#include "utils/genotype_reader.hpp"
#include "exceptions/program_error.hpp"

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

FacetWrapper FacetFactory::make(const std::string& name, const std::vector<VcfRecord>& records) const
{
    if (facet_makers_.count(name) == 1) {
        return facet_makers_.at(name)(records);
    } else {
        throw UnknownFacet {name};
    }
}

void FacetFactory::setup_facet_makers()
{
    facet_makers_["ReadAssignments"] = [this] (const std::vector<VcfRecord>& records) -> FacetWrapper
    {
        auto samples = read_pipe_.source().samples();
        auto genotypes = extract_genotypes(records, samples, reference_);
        ReadMap reads {};
        if (!records.empty()) {
            reads = read_pipe_.fetch_reads(encompassing_region(records));
        }
        return FacetWrapper {std::make_unique<ReadAssignments>(genotypes, reads)};
    };
}

} // namespace csr
} // namespace octopus
