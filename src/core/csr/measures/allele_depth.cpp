// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_depth.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/tools/read_assigner.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/genotype_reader.hpp"
#include "../facets/samples.hpp"
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
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments"));
    Array<Optional<Array<ValueType>>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample = samples[s];
        if (is_evaluable(call, sample)) {
            const auto alleles = get_called_alt_alleles(call, sample);
            const auto allele_support = compute_allele_support(alleles, assignments, sample);
            result[s] = Array<ValueType>(alleles.size());
            std::transform(std::cbegin(alleles), std::cend(alleles), std::begin(*result[s]),
                           [&] (const auto& allele) { return allele_support.at(allele).size(); });
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
    return {"Samples", "ReadAssignments"};
}

} // namespace csr
} // namespace octopus
