// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "contamination_sum.hpp"

#include <algorithm>
#include <iterator>
#include <set>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

#include "core/types/allele.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/tools/read_assigner.hpp"
#include "io/variant/vcf_record.hpp"
#include "contamination.hpp"

namespace octopus { namespace csr {

const std::string ContaminationSum::name_ = "CFS";

std::unique_ptr<Measure> ContaminationSum::do_clone() const
{
    return std::make_unique<ContaminationSum>(*this);
}

Measure::ValueType ContaminationSum::get_value_type() const
{
    return double {};
}

Measure::ResultType ContaminationSum::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    const auto contaminations = boost::get<Array<Array<Optional<ValueType>>>>(Contamination{}.evaluate(call, facets));
    const auto num_alleles = 1 + call.num_alt();
    Array<Optional<ValueType>> result(num_alleles);
    for (std::size_t a {0}; a < num_alleles; ++a) {
        for (std::size_t s {0}; s < contaminations.size(); ++s) {
            if (contaminations[s][a]) {
                if (result[a]) {
                    result[a] = boost::get<double>(*result[a]) + boost::get<double>(*contaminations[s][a]);
                } else {
                    result[a] = boost::get<double>(*contaminations[s][a]);
                }
            }
        }
    }
    return result;
}

Measure::ResultCardinality ContaminationSum::do_cardinality() const noexcept
{
    return ResultCardinality::alleles;
}

const std::string& ContaminationSum::do_name() const
{
    return name_;
}

std::string ContaminationSum::do_describe() const
{
    return "Sum of fraction of reads that support a haplotype called in another sample containing each the uncalled allele";
}

std::vector<std::string> ContaminationSum::do_requirements() const
{
    return Contamination{}.requirements();
}

boost::optional<Measure::Aggregator> ContaminationSum::do_aggregator() const noexcept
{
    return Measure::Aggregator::sum;
}

} // namespace csr
} // namespace octopus
