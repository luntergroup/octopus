// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "measure.hpp"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <type_traits>

#include <boost/lexical_cast.hpp>

#include "io/variant/vcf_spec.hpp"
#include "utils/append.hpp"

namespace octopus { namespace csr {

void Measure::do_set_parameters(std::vector<std::string> params)
{
    if (!params.empty()) {
        throw BadMeasureParameters {this->name()};
    }
}

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

struct MeasureResultTypeVisitor : boost::static_visitor<std::string>
{
    auto operator()(bool) const { return vcfspec::header::meta::type::flag; }
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
    auto operator()(T) const { return vcfspec::header::meta::type::integer; }
    template <typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
    auto operator()(T) const { return vcfspec::header::meta::type::floating; }
    template <typename T> auto operator()(boost::optional<T>) const { return (*this)(T{}); }
    template <typename T> auto operator()(std::vector<T>) const { return (*this)(T{}); }
    auto operator()(boost::any) const { return vcfspec::header::meta::type::string; }
};

std::string get_vcf_typename(const Measure::ResultType& value)
{
    MeasureResultTypeVisitor vis {};
    return boost::apply_visitor(vis, value);
}

void Measure::annotate(VcfHeader::Builder& header) const
{
    if (!is_required_vcf_field()) {
        const auto vcf_typename = get_vcf_typename(this->get_default_result());
        using namespace vcfspec::header::meta::number;
        if (this->cardinality() == Measure::ResultCardinality::num_samples) {
            header.add_format(this->name(), one, vcf_typename, this->describe());
        } else if (this->cardinality() == Measure::ResultCardinality::num_alleles) {
            header.add_info(this->name(), per_allele, vcf_typename, this->describe());
        } else {
            header.add_info(this->name(), one, vcf_typename, this->describe());
        }
    }
}

struct VectorIndexGetterVisitor : public boost::static_visitor<Measure::ResultType>
{
    VectorIndexGetterVisitor(std::size_t idx) : idx_ {idx} {}
    template <typename T> T operator()(const std::vector<T>& value) const noexcept { return value[idx_]; }
    template <typename T> boost::optional<T> operator()(const boost::optional<std::vector<T>>& value) const noexcept
    {
        if (value) {
            return (*value)[idx_];
        } else {
            return boost::none;
        }
    }
    template <typename T> T operator()(const T& value) const noexcept { return value; }
private:
    std::size_t idx_;
};

void Measure::annotate(VcfRecord::Builder& record, const ResultType& value, const VcfHeader& header) const
{
    if (!is_required_vcf_field()) {
        if (this->cardinality() == Measure::ResultCardinality::num_samples) {
            record.add_format(this->name());
            const auto samples = header.samples();
            for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
                const auto sample_value = boost::apply_visitor(VectorIndexGetterVisitor {sample_idx}, value);
                record.set_format(samples[sample_idx], this->name(), this->serialise(sample_value));
            }
        } else {
            record.set_info(this->name(), this->serialise(value));
        }
    }
}

// non-member methods

std::string long_name(const Measure& measure)
{
    auto result = measure.name();
    const auto params = measure.parameters();
    if (!params.empty()) {
        result += '[' + utils::join(params, ',') + ']';
    }
    return result;
}

std::string long_name(const MeasureWrapper& measure)
{
    return long_name(*measure.base());
}

struct IsMissingMeasureVisitor : public boost::static_visitor<bool>
{
    template <typename T> bool operator()(const boost::optional<T>& value) const noexcept { return !value; }
    template <typename T> bool operator()(const T& value) const noexcept { return false; }
};

bool is_missing(const Measure::ResultType& value) noexcept
{
    return boost::apply_visitor(IsMissingMeasureVisitor {}, value);
}

std::vector<std::string> get_all_requirements(const std::vector<MeasureWrapper>& measures)
{
    std::vector<std::string> result {};
    result.reserve(3 * measures.size()); // Just a guess
    for (const auto& measure : measures) {
        utils::append(measure.requirements(), result);
    }
    std::sort(std::begin(result), std::end(result));
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

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
