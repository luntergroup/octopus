// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "measure.hpp"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <cassert>

#include <boost/lexical_cast.hpp>

#include "io/variant/vcf_spec.hpp"

namespace octopus { namespace csr {

struct MeasureSerialiseVisitor : boost::static_visitor<>
{
    std::string str;
    
    void operator()(double value)
    {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(3) << value;
        str = ss.str();
    }
    void operator()(boost::any value)
    {
        str = ".";
    }
    template <typename T>
    void operator()(const T& value)
    {
        str = boost::lexical_cast<std::string>(value);
    }
    template <typename T>
    void operator()(const boost::optional<T>& value)
    {
        if (value) {
            (*this)(*value);
        } else {
            str = ".";
        }
    }
    template <typename T>
    void operator()(const std::vector<T>& values)
    {
        if (values.empty()) {
            str = ".";
        } else {
            auto tmp_str = std::move(str);
            std::for_each(std::cbegin(values), std::prev(std::cend(values)), [&] (const auto& value) {
                (*this)(value);
                tmp_str += str;
                tmp_str += ',';
            });
            (*this)(values.back());
            tmp_str += str;
            str = std::move(tmp_str);
        }
    }
};

std::string Measure::do_serialise(const ResultType& value) const
{
    MeasureSerialiseVisitor vis {};
    boost::apply_visitor(vis, value);
    return vis.str;
}

void Measure::annotate(VcfHeader::Builder& header) const
{
    if (!is_required_vcf_field()) {
        std::string number;
        using namespace vcfspec::header::meta::number;
        switch (this->cardinality()) {
            case Measure::ResultCardinality::num_samples: number = unknown; break;
            case Measure::ResultCardinality::num_alleles: number = per_allele; break;
            case Measure::ResultCardinality::one: number = "1"; break;
        }
        header.add_info(this->name(), number, "String", this->describe());
    }
}

void Measure::annotate(VcfRecord::Builder& record, const ResultType& value) const
{
    if (!is_required_vcf_field()) {
        record.set_info(this->name(), this->serialise(value));
    }
}

// non-member methods

struct IsMissingMeasureVisitor : public boost::static_visitor<bool>
{
    template <typename T> bool operator()(const boost::optional<T>& value) const noexcept { return !value; }
    template <typename T> bool operator()(const T& value) const noexcept { return false; }
};

bool is_missing(const Measure::ResultType& value) noexcept
{
    return boost::apply_visitor(IsMissingMeasureVisitor {}, value);
}

struct VectorIndexGetterVisitor : public boost::static_visitor<Measure::ResultType>
{
    VectorIndexGetterVisitor(std::size_t idx) : idx_ {idx} {}
    template <typename T> T operator()(const std::vector<T>& value) const noexcept { return value[idx_]; }
    template <typename T> T operator()(const T& value) const noexcept { return value; }
private:
    std::size_t idx_;
};

Measure::ResultType get_sample_value(const Measure::ResultType& value, const MeasureWrapper& measure, const std::size_t sample_idx)
{
    if (measure.cardinality() == Measure::ResultCardinality::num_samples) {
        return boost::apply_visitor(VectorIndexGetterVisitor {sample_idx}, value);
    } else {
        return value;
    }
}

std::vector<Measure::ResultType>
get_sample_values(const std::vector<Measure::ResultType>& values, const std::vector<MeasureWrapper>& measures, std::size_t sample_idx)
{
    assert(values.size() == measures.size());
    std::vector<Measure::ResultType> result(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::cbegin(measures), std::begin(result),
                   [&] (const auto& value, const auto& measure) { return get_sample_value(value, measure, sample_idx); });
    return result;
}

} // namespace csr
} // namespace octopus
