// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "allele_mapping_quality.hpp"

#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <numeric>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include "basics/aligned_read.hpp"
#include "core/types/allele.hpp"
#include "io/variant/vcf_record.hpp"
#include "utils/maths.hpp"
#include "../facets/samples.hpp"
#include "../facets/alleles.hpp"
#include "../facets/read_assignments.hpp"

namespace octopus { namespace csr {

const std::string AlleleMappingQuality::name_ = "AMQ";

std::unique_ptr<Measure> AlleleMappingQuality::do_clone() const
{
    return std::make_unique<AlleleMappingQuality>(*this);
}

Measure::ValueType AlleleMappingQuality::get_value_type() const
{
    return int {};
}

namespace {

int calculate_median_mapping_quality(const ReadRefSupportSet& reads)
{
    std::vector<AlignedRead::MappingQuality> mapping_qualities(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(mapping_qualities),
                   [](const AlignedRead& read) { return read.mapping_quality(); });
    return maths::median(mapping_qualities);
}

} // namespace

Measure::ResultType AlleleMappingQuality::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
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
                    result[s][a] = calculate_median_mapping_quality(support_set_itr->second);
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality AlleleMappingQuality::do_cardinality() const noexcept
{
    return ResultCardinality::samples_and_alleles;
}

const std::string& AlleleMappingQuality::do_name() const
{
    return name_;
}

std::string AlleleMappingQuality::do_describe() const
{
    return "Median mapping quality of reads assigned to each allele";
}

std::vector<std::string> AlleleMappingQuality::do_requirements() const
{
    return {"Samples", "Alleles", "ReadAssignments"};
}

boost::optional<Measure::Aggregator> AlleleMappingQuality::do_aggregator() const noexcept
{
    return Measure::Aggregator::min_tail;
}

} // namespace csr
} // namespace octopus
