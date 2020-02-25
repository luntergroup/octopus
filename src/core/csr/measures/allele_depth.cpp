// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_depth.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

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

Measure::ResultType AlleleDepth::get_default_result() const
{
    return std::vector<boost::optional<int>> {};
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

template <typename T>
void pop_front(std::vector<T>& v)
{
    v.erase(std::cbegin(v));
}

auto get_support_counts(const std::vector<Allele>& alleles, const AlleleSupportMap& support)
{
    std::vector<int> result(alleles.size());
    std::transform(std::cbegin(alleles), std::cend(alleles), std::begin(result),
                   [&] (const auto& allele) { return support.at(allele).size(); });
    return result;
}

auto min_support_count(const std::vector<Allele>& alleles, const AlleleSupportMap& support)
{
    const auto counts = get_support_counts(alleles, support);
    return *std::min_element(std::cbegin(counts), std::cend(counts));
}

} // namespace

Measure::ResultType AlleleDepth::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    std::vector<boost::optional<int>> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        boost::optional<int> sample_result {};
        if (is_evaluable(call, sample)) {
            std::vector<Allele> alleles; bool has_ref;
            std::tie(alleles, has_ref) = get_called_alleles(call, sample);
            if (has_ref) pop_front(alleles); // ref always first
            sample_result = min_support_count(alleles, assignments.at(sample));
        }
        result.push_back(sample_result);
    }
    return result;
}

Measure::ResultCardinality AlleleDepth::do_cardinality() const noexcept
{
    return ResultCardinality::samples;
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
