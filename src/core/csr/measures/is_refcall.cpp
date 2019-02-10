// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "is_refcall.hpp"

#include "io/variant/vcf_record.hpp"
#include "config/octopus_vcf.hpp"
#include "../facets/samples.hpp"

namespace octopus { namespace csr {

const std::string IsRefcall::name_ = "REFCALL";

IsRefcall::IsRefcall(bool report_sample_status) : report_sample_status_ {report_sample_status} {}

std::unique_ptr<Measure> IsRefcall::do_clone() const
{
    return std::make_unique<IsRefcall>(*this);
}

namespace {

bool is_refcall(const VcfRecord& record)
{
    return record.alt().size() == 1 && record.alt()[0] == vcf::spec::allele::nonref;
}

} // namespace

Measure::ResultType IsRefcall::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (report_sample_status_) {
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        std::vector<bool> result(samples.size());
        std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                       [&call] (const auto& sample) { return call.is_homozygous_ref(sample); });
        return result;
    } else {
        return is_refcall(call);
    }
}

Measure::ResultCardinality IsRefcall::do_cardinality() const noexcept
{
    if (report_sample_status_) {
        return ResultCardinality::num_samples;
    } else {
        return ResultCardinality::one;
    }
}

const std::string& IsRefcall::do_name() const
{
    return name_;
}

std::string IsRefcall::do_describe() const
{
    if (report_sample_status_) {
        return "REFCALL status of each sample";
    } else {
        return "Is the call marked REFCALL";
    }
}

std::vector<std::string> IsRefcall::do_requirements() const
{
    if (report_sample_status_) {
        return {"Samples"};
    } else {
        return {};
    }
}

bool IsRefcall::is_equal(const Measure& other) const noexcept
{
    return report_sample_status_ == static_cast<const IsRefcall&>(other).report_sample_status_;
}

} // namespace csr
} // namespace octopus
