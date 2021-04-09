// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_frequency_bias.hpp"

#include <algorithm>
#include <iterator>
#include <numeric>
#include <cmath>

#include <boost/variant.hpp>

#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/concat.hpp"
#include "../facets/samples.hpp"
#include "alt_allele_count.hpp"
#include "allele_frequency.hpp"

namespace octopus { namespace csr {

const std::string AlleleFrequencyBias::name_ = "AFB";

std::unique_ptr<Measure> AlleleFrequencyBias::do_clone() const
{
    return std::make_unique<AlleleFrequencyBias>(*this);
}

Measure::ValueType AlleleFrequencyBias::get_value_type() const
{
    return double {};
}

namespace {

auto sum_alt_allele_counts(const Measure::Array<Measure::ValueType>& counts)
{
    return std::accumulate(std::cbegin(counts), std::cend(counts), 0u, [] (auto total, auto count) { return total + boost::get<int>(count); });
}

} // namespace

Measure::ResultType AlleleFrequencyBias::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto alt_allele_counts = boost::get<Array<Array<ValueType>>>(AltAlleleCount{}.evaluate(call, facets));
    const auto allele_frequencies = boost::get<Array<Array<Optional<ValueType>>>>(AlleleFrequency{}.evaluate(call, facets));
    const auto num_alleles = 1 + call.alt().size();
    const auto num_samples = samples.size();
    assert(allele_frequencies.size() == num_samples);
    assert(alt_allele_counts.size() == num_samples);
    Array<Array<Optional<ValueType>>> result(num_samples, Array<Optional<ValueType>>(num_alleles));
    for (std::size_t s {0}; s < num_samples; ++s) {
        assert(alt_allele_counts[s].size() == num_alleles - 1);
        assert(allele_frequencies[s].size() == num_alleles);
        const auto ploidy = call.genotype(samples[s]).size();
        const auto non_ref_count = sum_alt_allele_counts(alt_allele_counts[s]);
        assert(ploidy >= non_ref_count);
        const auto ref_count = ploidy - non_ref_count;
        for (std::size_t a {0}; a < num_alleles; ++a) {
            if (allele_frequencies[s][a]) {
                const auto expected_frequency = static_cast<double>(a == 0 ? ref_count : boost::get<int>(alt_allele_counts[s][a - 1])) / ploidy;
                const auto emperical_frequency = boost::get<double>(*allele_frequencies[s][a]);
                result[s][a] = std::abs(expected_frequency - emperical_frequency);
            }
        }
    }
    return result;
}

Measure::ResultCardinality AlleleFrequencyBias::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
}

const std::string& AlleleFrequencyBias::do_name() const
{
    return name_;
}

std::string AlleleFrequencyBias::do_describe() const
{
    return "Absolute difference between empirical allele frequency (AF) and expected allele frequency given genotype";
}

std::vector<std::string> AlleleFrequencyBias::do_requirements() const
{
    auto result = concat(AltAlleleCount{}.requirements(), AlleleFrequency{}.requirements());
    result.push_back("Samples");
    return result;
}

boost::optional<Measure::Aggregator> AlleleFrequencyBias::do_aggregator() const noexcept
{
    return Measure::Aggregator::min_tail;
}

} // namespace csr
} // namespace octopus
