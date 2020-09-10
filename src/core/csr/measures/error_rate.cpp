// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "error_rate.hpp"

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
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string ErrorRate::name_ = "ER";

std::unique_ptr<Measure> ErrorRate::do_clone() const
{
    return std::make_unique<ErrorRate>(*this);
}

Measure::ValueType ErrorRate::get_value_type() const
{
    return double {};
}

namespace {

boost::optional<double> 
calculate_error_rate(const ReadRefSupportSet& reads)
{
    std::size_t error_bases {0}, total_bases {0};
    for (const AlignedRead& read : reads) {
        error_bases += sum_non_matches(read.cigar());
        total_bases += sequence_size(read);
    }
    if (total_bases > 0) {
        return static_cast<double>(error_bases) / total_bases;
    } else {
        return boost::none;
    }
}

} // namespace

Measure::ResultType ErrorRate::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
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
                    result[s][a] = calculate_error_rate(support_set_itr->second);
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality ErrorRate::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
}

const std::string& ErrorRate::do_name() const
{
    return name_;
}

std::string ErrorRate::do_describe() const
{
    return "Error rate in supporting reads";
}

std::vector<std::string> ErrorRate::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

boost::optional<Measure::Aggregator> ErrorRate::do_aggregator() const noexcept
{
    return Measure::Aggregator::max;
}
    
} // namespace csr
} // namespace octopus
