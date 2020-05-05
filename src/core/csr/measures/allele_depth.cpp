// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_depth.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string AlleleDepth::name_ = "AD";

std::unique_ptr<Measure> AlleleDepth::do_clone() const
{
    return std::make_unique<AlleleDepth>(*this);
}

Measure::ValueType AlleleDepth::get_value_type() const
{
    return std::size_t {};
}

namespace {

bool is_canonical(const VcfRecord::NucleotideSequence& allele) noexcept
{
    const static VcfRecord::NucleotideSequence deleted_allele {vcfspec::deletedBase};
    return !(allele == vcfspec::missingValue || allele == deleted_allele);
}

bool has_called_alt_allele(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    if (!call.has_genotypes()) return true;
    const auto& genotype = get_genotype(call, sample);
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                       [&] (const auto& allele) { return allele != call.ref() && is_canonical(allele); });
}

bool is_evaluable(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    return has_called_alt_allele(call, sample);
}

} // namespace

Measure::ResultType AlleleDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    Array<Optional<Array<ValueType>>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample = samples[s];
        if (is_evaluable(call, sample)) {
            const auto sample_alleles = get_all(alleles, call, sample);
            result[s] = Array<ValueType>(sample_alleles.size());
            const auto& support = assignments.at(sample);
            std::transform(std::cbegin(sample_alleles), std::cend(sample_alleles), std::begin(*result[s]),
                           [&] (const auto& allele) { return support.at(allele).size(); });
        }
    }
    return result;
}

Measure::ResultCardinality AlleleDepth::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alt_alleles;
}

const std::string& AlleleDepth::do_name() const
{
    return name_;
}

std::string AlleleDepth::do_describe() const
{
    return "Minor empirical alt allele depth";
}

std::vector<std::string> AlleleDepth::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
