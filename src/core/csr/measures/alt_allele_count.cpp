// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "alt_allele_count.hpp"

#include <iterator>
#include <algorithm>

#include <boost/variant.hpp>

#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "../facets/samples.hpp"

namespace octopus { namespace csr {

const std::string AltAlleleCount::name_ = "AC";

std::unique_ptr<Measure> AltAlleleCount::do_clone() const
{
    return std::make_unique<AltAlleleCount>(*this);
}

namespace {

int count_non_ref_alleles(const VcfRecord& call, const VcfRecord::SampleName& sample)
{
    const auto genotype = get_genotype(call, sample);
    return genotype.size() - std::count(std::cbegin(genotype), std::cend(genotype), call.ref());
}

} // namespace

Measure::ResultType AltAlleleCount::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    std::vector<int> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        result.push_back(count_non_ref_alleles(call, sample));
    }
    return result;
}

Measure::ResultCardinality AltAlleleCount::do_cardinality() const noexcept
{
    return ResultCardinality::num_samples;
}

const std::string& AltAlleleCount::do_name() const
{
    return name_;
}

std::string AltAlleleCount::do_describe() const
{
    return "Number of non-reference alleles called for each sample";
}

std::vector<std::string> AltAlleleCount::do_requirements() const
{
    return {"Samples"};
}

} // namespace csr
} // namespace octopus
