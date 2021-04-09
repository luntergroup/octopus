// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_rate_stdev.hpp"

#include <numeric>
#include <algorithm>
#include <iterator>
#include <map>
#include <cstddef>

#include <boost/variant.hpp>

#include "basics/cigar_string.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/haplotype.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ErrorRateStdev::name_ = "ERS";

std::unique_ptr<Measure> ErrorRateStdev::do_clone() const
{
    return std::make_unique<ErrorRateStdev>(*this);
}

Measure::ValueType ErrorRateStdev::get_value_type() const
{
    return double {};
}

namespace {

auto calculate_error_rate(const AlignedRead& read) noexcept
{
    return static_cast<double>(sum_non_matches(read.cigar())) / sequence_size(read);
}

boost::optional<double>
calculate_error_rate_stdev(const ReadRefSupportSet& reads)
{
    std::vector<double> error_rates(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(error_rates), calculate_error_rate);
    if (!error_rates.empty()) {
        return maths::stdev(error_rates);
    } else {
        return boost::none;
    }
}

} // namespace

Measure::ResultType ErrorRateStdev::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto& samples = get_value<Samples>(facets.at("Samples"));
    const auto& alleles = get_value<Alleles>(facets.at("Alleles"));
    const auto& assignments = get_value<ReadAssignments>(facets.at("ReadAssignments")).alleles;
    const auto num_alleles = 1 + call.num_alt();
    Array<Array<Optional<ValueType>>> result(samples.size(), Array<Optional<ValueType>>(num_alleles));
    for (std::size_t s {0}; s < samples.size(); ++s) {
        const auto& sample_alleles = get(alleles, call, samples[s]);
        assert(sample_alleles.size() == num_alleles);
        const auto& support = assignments.at(samples[s]);
        for (std::size_t a {0}; a < num_alleles; ++a) {
            if (sample_alleles[a]) {
                const auto support_set_itr = support.find(*sample_alleles[a]);
                if (support_set_itr != std::cend(support)) {
                    result[s][a] = calculate_error_rate_stdev(support_set_itr->second);
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality ErrorRateStdev::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
}

const std::string& ErrorRateStdev::do_name() const
{
    return name_;
}

std::string ErrorRateStdev::do_describe() const
{
    return "Error rate standard deviation in supporting reads";
}

std::vector<std::string> ErrorRateStdev::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

boost::optional<Measure::Aggregator> ErrorRateStdev::do_aggregator() const noexcept
{
    return Measure::Aggregator::max;
}

} // namespace csr
} // namespace octopus
