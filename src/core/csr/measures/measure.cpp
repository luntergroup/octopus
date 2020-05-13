// Copyright (c) 2015-2020 Daniel Cooke
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
#include "utils/string_utils.hpp"

namespace octopus { namespace csr {

void Measure::do_set_parameters(std::vector<std::string> params)
{
    if (!params.empty()) {
        throw BadMeasureParameters {this->name()};
    }
}

struct MeasureSerialiseVisitor : boost::static_visitor<>
{
    std::ostringstream ss;
    
    void operator()(const Measure::ValueType& value)
    {
        boost::apply_visitor(*this, value);
    }
    void operator()(double value)
    {
        ss << std::fixed << std::setprecision(3) << value;
    }
    template <typename T>
    void operator()(const T& value)
    {
        ss << boost::lexical_cast<std::string>(value);
    }
    template <typename T>
    void operator()(const Measure::Optional<T>& value)
    {
        if (value) {
            (*this)(*value);
        } else {
            ss << vcfspec::missingValue;
        }
    }
    template <typename T>
    void operator()(const Measure::Array<T>& values)
    {
        if (values.empty()) {
            ss << vcfspec::missingValue;
        } else {
            std::for_each(std::cbegin(values), std::prev(std::cend(values)), [&] (const auto& value) {
                (*this)(value);
                ss << vcfspec::format::valueSeperator;
            });
            (*this)(values.back());
        }
    }
};

std::string Measure::do_serialise(const ResultType& value) const
{
    MeasureSerialiseVisitor vis {};
    boost::apply_visitor(vis, value);
    return vis.ss.str();
}

struct MeasureValueTypeVisitor : boost::static_visitor<std::string>
{
    auto operator()(bool) const { return vcfspec::header::meta::type::flag; }
    template <typename T, std::enable_if_t<std::is_integral<T>::value, int> = 0>
    auto operator()(T) const { return vcfspec::header::meta::type::integer; }
    template <typename T, std::enable_if_t<std::is_floating_point<T>::value, int> = 0>
    auto operator()(T) const { return vcfspec::header::meta::type::floating; }
    template <typename T> auto operator()(boost::optional<T>) const { return (*this)(T{}); }
    template <typename T> auto operator()(std::vector<T>) const { return (*this)(T{}); }
};

std::string get_vcf_typename(const Measure::ValueType& value)
{
    MeasureValueTypeVisitor vis {};
    return boost::apply_visitor(vis, value);
}

bool is_per_sample(const Measure::ResultCardinality cardinality) noexcept
{
    using RC = Measure::ResultCardinality;
    return cardinality == RC::samples
        || cardinality == RC::samples_and_alleles
        || cardinality == RC::samples_and_alt_alleles;
}

bool is_one_value_annotation(const Measure::ResultCardinality cardinality) noexcept
{
    using RC = Measure::ResultCardinality;
    return cardinality == RC::one || cardinality == RC::samples;
}

auto get_vcf_number(const Measure::ResultCardinality cardinality, const bool aggregate_alleles)
{
    using RC = Measure::ResultCardinality;
    using namespace vcfspec::header::meta::number;
    if (cardinality == RC::one || cardinality == RC::samples || aggregate_alleles) {
        return one;
    } else if (cardinality == RC::alleles || cardinality == RC::samples_and_alleles) {
        return per_allele;
    } else {
        return per_alt_allele;
    }
}

void Measure::annotate(VcfHeader::Builder& header, const bool aggregate_alleles) const
{
    if (!is_required_vcf_field()) {
        const auto vcf_typename = get_vcf_typename(this->get_value_type());
        const auto vcf_number = get_vcf_number(this->cardinality(), aggregate_alleles);
        if (is_per_sample(this->cardinality())) {
            header.add_format(this->name(), vcf_number, vcf_typename, this->describe());
        } else {
            header.add_info(this->name(), vcf_number, vcf_typename, this->describe());
        }
    }
}

struct VectorIndexGetterVisitor : public boost::static_visitor<Measure::ResultType>
{
    VectorIndexGetterVisitor(std::size_t idx) : idx_ {idx} {}
    template <typename T> Measure::ResultType operator()(const T& value) const 
    {
        return value;
    }
    template <typename T> Measure::ResultType operator()(const Measure::Array<T>& value) const
    {
         return value[idx_];
    }
    template <typename T> Measure::ResultType operator()(const Measure::Optional<Measure::Array<T>>& value) const
    {
        Measure::Optional<T> result {};
        if (value) {
            result = (*value)[idx_];
        }
        return result;
    }
    template <typename T> Measure::ResultType operator()(const Measure::Optional<Measure::Array<Measure::Optional<T>>>& value) const
    {
        Measure::Optional<T> result {};
        if (value) {
            result = (*value)[idx_];
        }
        return result;
    }
private:
    std::size_t idx_;
};

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

namespace {

struct PlusVisitor : public boost::static_visitor<Measure::ValueType>
{
    PlusVisitor(Measure::ValueType other) : other_ {other} {}
    template <typename T> auto operator()(const T value) const noexcept
    {
        return boost::get<T>(other_) + value;
    }
private:
    Measure::ValueType other_;
};

auto sum(const Measure::Array<Measure::ValueType>& values)
{
    return std::accumulate(std::next(std::cbegin(values)), std::cend(values), values.front(),
                           [] (auto total, const auto& value) { return boost::apply_visitor(PlusVisitor {total}, value); });
}

auto min(const Measure::Array<Measure::ValueType>& values)
{
    return *std::min_element(std::cbegin(values), std::cend(values));
}

auto max(const Measure::Array<Measure::ValueType>& values)
{
    return *std::max_element(std::cbegin(values), std::cend(values));
}

struct DivideVisitor : public boost::static_visitor<Measure::ValueType>
{
    DivideVisitor(std::size_t n) : n_ {n} {}
    template <typename T> auto operator()(const T value) const noexcept { return value / n_; }
private:
    std::size_t n_;
};

auto mean(const Measure::Array<Measure::ValueType>& values)
{
    return boost::apply_visitor(DivideVisitor {values.size()}, sum(values));
}

Measure::Optional<Measure::ValueType> 
aggregate(const Measure::Array<Measure::ValueType>& values,
          const Measure::Aggregator aggregator)
{
    if (values.empty()) {
        return boost::none;
    } else {
        switch (aggregator) {
            case Measure::Aggregator::min: return min(values);
            case Measure::Aggregator::max: return max(values);
            case Measure::Aggregator::sum: return sum(values);
            case Measure::Aggregator::mean: return mean(values);
            default: return boost::none;
        }
    }
}

} // namespace

struct AggregatorVisitor : public boost::static_visitor<Measure::Optional<Measure::ValueType>>
{
    AggregatorVisitor(Measure::Aggregator aggregator) : aggregator_ {aggregator} {}
    Measure::Optional<Measure::ValueType> operator()(const Measure::ValueType& value) const 
    {
        return value;
    }
    Measure::Optional<Measure::ValueType> operator()(const Measure::Optional<Measure::ValueType>& value) const 
    {
        return value;
    }
    Measure::Optional<Measure::ValueType> operator()(const Measure::Array<Measure::ValueType>& values) const
    {
        return aggregate(values, aggregator_);
    }
    Measure::Optional<Measure::ValueType> operator()(const Measure::Array<Measure::Optional<Measure::ValueType>>& values) const
    {
        Measure::Array<Measure::ValueType> valid_values {};
        valid_values.reserve(values.size());
        for (const auto& value : values) {
            if (value) valid_values.push_back(*value);
        }
        return (*this)(valid_values);
    }
    template <typename T> 
    Measure::Optional<Measure::ValueType> operator()(const Measure::Optional<Measure::Array<T>>& values) const
    {
        if (values) {
            return (*this)(*values);
        } else {
            return boost::none;
        }
    }
    template <typename T> 
    Measure::Optional<Measure::ValueType> operator()(const T& values) const
    {
        throw std::runtime_error {"Bad AggregatorVisitor"};
    }
private:
    Measure::Aggregator aggregator_;
};

void Measure::annotate(VcfRecord::Builder& record, const ResultType& value, const VcfHeader& header, const bool aggregate_alleles) const
{
    if (!is_required_vcf_field()) {
        if (is_per_sample(this->cardinality())) {
            record.add_format(this->name());
            const auto samples = header.samples();
            for (std::size_t sample_idx {0}; sample_idx < samples.size(); ++sample_idx) {
                auto sample_value = boost::apply_visitor(VectorIndexGetterVisitor {sample_idx}, value);
                bool split_alleles_required {!is_one_value_annotation(this->cardinality())};
                if (aggregate_alleles && this->aggregator()) {
                    sample_value = boost::apply_visitor(AggregatorVisitor {*this->aggregator()}, sample_value);
                    split_alleles_required = false;
                }
                if (!split_alleles_required) {
                    record.set_format(samples[sample_idx], this->name(), this->serialise(sample_value));
                } else {
                    record.set_format(samples[sample_idx], this->name(), utils::split(this->serialise(sample_value), vcfspec::format::valueSeperator));
                }
            }
        } else {
            if (aggregate_alleles && this->aggregator()) {
                record.set_info(this->name(), this->serialise(boost::apply_visitor(AggregatorVisitor {*this->aggregator()}, value)));
            } else {
                record.set_info(this->name(), this->serialise(value));
            }
        }
    }
}

// non-member methods

Measure::ResultType get_sample_value(const Measure::ResultType& value, const MeasureWrapper& measure, const std::size_t sample_idx, const bool aggregate)
{
    Measure::ResultType result {};
    if (is_per_sample(measure.cardinality())) {
        result = boost::apply_visitor(VectorIndexGetterVisitor {sample_idx}, value);
    } else {
        result = value;
    }
    if (aggregate && measure.aggregator()) {
        result = boost::apply_visitor(AggregatorVisitor {*measure.aggregator()}, result);
    }
    return result;
}

std::vector<Measure::ResultType>
get_sample_values(const std::vector<Measure::ResultType>& values, const std::vector<MeasureWrapper>& measures, std::size_t sample_idx, const bool aggregate)
{
    assert(values.size() == measures.size());
    std::vector<Measure::ResultType> result(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::cbegin(measures), std::begin(result),
                   [&] (const auto& value, const auto& measure) { return get_sample_value(value, measure, sample_idx, aggregate); });
    return result;
}

} // namespace csr
} // namespace octopus
