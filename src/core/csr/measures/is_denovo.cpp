// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "is_denovo.hpp"

#include "basics/pedigree.hpp"
#include "io/variant/vcf_record.hpp"
#include "config/octopus_vcf.hpp"
#include "../facets/samples.hpp"
#include "../facets/pedigree.hpp"

namespace octopus { namespace csr {

const std::string IsDenovo::name_ = "DENOVO";

IsDenovo::IsDenovo(bool report_sample_status) : report_sample_status_ {report_sample_status} {}

std::unique_ptr<Measure> IsDenovo::do_clone() const
{
    return std::make_unique<IsDenovo>(*this);
}

namespace {

bool is_denovo(const VcfRecord& call)
{
    return call.has_info(vcf::spec::info::denovo);
}

auto child_idx(const std::vector<SampleName>& samples, const octopus::Pedigree& pedigree)
{
    assert(samples.size() == 3);
    if (is_parent_of(samples[0], samples[1], pedigree)) {
        return 1;
    } else if (is_parent_of(samples[1], samples[0], pedigree)) {
        return 0;
    } else {
        return 2;
    }
}

} // namespace

Measure::ResultType IsDenovo::do_evaluate(const VcfRecord& call, const FacetMap& facets) const
{
    if (report_sample_status_) {
        const auto& samples = get_value<Samples>(facets.at("Samples"));
        std::vector<bool> result(samples.size(), false);
        if (is_denovo(call)) {
            const auto& pedigree = get_value<Pedigree>(facets.at("Pedigree"));
            result[child_idx(samples, pedigree)] = true;
        }
        return result;
    } else {
        return is_denovo(call);
    }
}

Measure::ResultCardinality IsDenovo::do_cardinality() const noexcept
{
    if (report_sample_status_) {
        return ResultCardinality::num_samples;
    } else {
        return ResultCardinality::one;
    }
}

const std::string& IsDenovo::do_name() const
{
    return name_;
}

std::string IsDenovo::do_describe() const
{
    if (report_sample_status_) {
        return "DENOVO status of each sample";
    } else {
        return "Is the call marked DENOVO";
    }
}

std::vector<std::string> IsDenovo::do_requirements() const
{
    if (report_sample_status_) {
        return {"Samples", "Pedigree"};
    } else {
        return {};
    }
}

bool IsDenovo::is_equal(const Measure& other) const noexcept
{
    return report_sample_status_ == static_cast<const IsDenovo&>(other).report_sample_status_;
}

} // namespace csr
} // namespace octopus
