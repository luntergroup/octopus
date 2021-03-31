// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "mean_likelihood.hpp"

#include <algorithm>
#include <iterator>

#include <boost/variant.hpp>

#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "io/variant/vcf_spec.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/genotypes.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string MeanLikelihood::name_ = "MHL";

std::unique_ptr<Measure> MeanLikelihood::do_clone() const
{
    return std::make_unique<MeanLikelihood>(*this);
}

Measure::ValueType MeanLikelihood::get_value_type() const
{
    return int {};
}

namespace {

const auto& get_genotype(const Facet::GenotypeMap& genotypes, const VcfRecord& call, const SampleName& sample)
{
    return overlap_range(genotypes.at(sample), call).front();
}

} // namespace

Measure::ResultType MeanLikelihood::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).haplotypes;
    const auto& genotypes = get_value<Genotypes>(facets.at("Genotypes"));
    Array<Array<Optional<ValueType>>> result(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample = samples[s];
        const auto& genotype = get_genotype(genotypes, call, sample);
        assert(genotype.ploidy() == call.ploidy(sample));
        result[s] = Array<Optional<ValueType>>(genotype.ploidy());
        for (unsigned p {0}; p < genotype.ploidy(); ++p) {
            const auto& assigned_reads = assignments.at(sample).assigned_wrt_haplotype.at(genotype[p]);
            const auto& assigned_likelihoods = assignments.at(sample).assigned_likelihoods.at(genotype[p]);
            assert(assigned_reads.size() == assigned_likelihoods.size());
            int phreds {0}, count {0};
            for (std::size_t n {0}; n < assigned_reads.size(); ++n) {
                if (overlaps(assigned_reads[n], call)) {
                    phreds += assigned_likelihoods[n] / -maths::constants::ln10Div10<>;
                    ++count;
                }
            }
            if (count > 0) result[s][p] = phreds / count;
        }
    }
    return result;
}

Measure::ResultCardinality MeanLikelihood::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_haplotypes;
}

const std::string& MeanLikelihood::do_name() const
{
    return name_;
}

std::string MeanLikelihood::do_describe() const
{
    return "Mean likelihood (Phreds) of reads overlapping the site assigned to each haplotype";
}

std::vector<std::string> MeanLikelihood::do_requirements() const
{
    return {"Samples", "Genotypes", "ReadAssignments"};
}

boost::optional<Measure::Aggregator> MeanLikelihood::do_aggregator() const noexcept
{
    return Measure::Aggregator::max;
}

} // namespace csr
} // namespace octopus
