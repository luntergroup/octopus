// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef threshold_filter_hpp
#define threshold_filter_hpp

#include <vector>
#include <string>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include "variant_call_filter.hpp"
#include "logging/progress_meter.hpp"
#include "../facets/facet_factory.hpp"

namespace octopus {

class VcfHeader;

namespace csr {

class ThresholdVariantCallFilter : public VariantCallFilter
{
public:
    using VariantCallFilter::OutputOptions;
    
    struct Threshold
    {
        virtual ~Threshold() = default;
        virtual std::unique_ptr<Threshold> clone() const = 0;
        virtual bool operator()(Measure::ResultType value) const noexcept = 0;
    };
    
    struct ThresholdWrapper
    {
        ThresholdWrapper(std::unique_ptr<Threshold> threshold) : threshold {std::move(threshold)} {}
        ThresholdWrapper(const ThresholdWrapper& other) : threshold {other.threshold->clone()} {}
        ThresholdWrapper& operator=(const ThresholdWrapper& other)
        {
            if (&other != this) threshold = other.threshold->clone();
            return *this;
        }
        bool operator()(Measure::ResultType value) const noexcept { return (*threshold)(value); }
        std::unique_ptr<Threshold> threshold;
    };
    
    struct Condition
    {
        MeasureWrapper measure;
        ThresholdWrapper threshold;
        std::string vcf_filter_key = ".";
    };
    
    ThresholdVariantCallFilter() = delete;
    
    ThresholdVariantCallFilter(FacetFactory facet_factory,
                               std::vector<Condition> conditions,
                               OutputOptions output_config,
                               boost::optional<ProgressMeter&> progress = boost::none);
    
    ThresholdVariantCallFilter(const ThresholdVariantCallFilter&)            = delete;
    ThresholdVariantCallFilter& operator=(const ThresholdVariantCallFilter&) = delete;
    ThresholdVariantCallFilter(ThresholdVariantCallFilter&&)                 = default;
    ThresholdVariantCallFilter& operator=(ThresholdVariantCallFilter&&)      = default;
    
    virtual ~ThresholdVariantCallFilter() = default;

private:
    std::vector<ThresholdWrapper> thresholds_;
    std::vector<std::string> vcf_filter_keys_;
    bool all_unique_filter_keys_;
    
    virtual void annotate(VcfHeader::Builder& header) const override;
    virtual Classification classify(const MeasureVector& measures) const override;
    
    bool passes_all_filters(const MeasureVector& measures) const;
    std::vector<std::string> get_failing_vcf_filter_keys(const MeasureVector& measures) const;
};

template <typename M, typename... Args>
decltype(auto) make_wrapped_threshold(Args&&... args)
{
    using Wrapper = ThresholdVariantCallFilter::ThresholdWrapper;
    return Wrapper {std::make_unique<M>(std::forward<Args>(args)...)};
}

struct LessThreshold : public ThresholdVariantCallFilter::Threshold
{
    explicit LessThreshold(double target) : visitor_ {target} {}
    std::unique_ptr<ThresholdVariantCallFilter::Threshold> clone() const
    {
        return std::make_unique<LessThreshold>(*this);
    }
    bool operator()(Measure::ResultType value) const noexcept
    {
        return boost::apply_visitor(visitor_, value);
    }
private:
    struct LessVisitor : public boost::static_visitor<bool>
    {
        explicit LessVisitor(double target) : target {target} {}
        bool operator()(double value) const noexcept { return value >= target; }
        bool operator()(boost::optional<double> value) const noexcept
        {
            return !value || (*this)(*value);
        }
        double target;
    };
    LessVisitor visitor_;
};

struct GreaterThreshold : public ThresholdVariantCallFilter::Threshold
{
    explicit GreaterThreshold(double target) : visitor_ {target} {}
    std::unique_ptr<ThresholdVariantCallFilter::Threshold> clone() const
    {
        return std::make_unique<GreaterThreshold>(*this);
    }
    bool operator()(Measure::ResultType value) const noexcept
    {
        return boost::apply_visitor(visitor_, value);
    }
private:
    struct GreaterVisitor : public boost::static_visitor<bool>
    {
        explicit GreaterVisitor(double target) : target {target} {}
        bool operator()(double value) const noexcept { return value <= target; }
        bool operator()(boost::optional<double> value) const noexcept
        {
            return !value || (*this)(*value);
        }
        double target;
    };
    GreaterVisitor visitor_;
};

struct BetweenThreshold : public ThresholdVariantCallFilter::Threshold
{
    explicit BetweenThreshold(double lower_bound, double upper_bound) : visitor_ {lower_bound, upper_bound} {}
    std::unique_ptr<ThresholdVariantCallFilter::Threshold> clone() const
    {
        return std::make_unique<BetweenThreshold>(*this);
    }
    bool operator()(Measure::ResultType value) const noexcept
    {
        return boost::apply_visitor(visitor_, value);
    }
private:
    struct BetweenVisitor : public boost::static_visitor<bool>
    {
        explicit BetweenVisitor(double lower_bound, double upper_bound)
        : lower_bound {lower_bound}, upper_bound {upper_bound} {}
        bool operator()(double value) const noexcept
        {
            return lower_bound <= value && value <= upper_bound;
        }
        bool operator()(boost::optional<double> value) const noexcept
        {
            return !value || (*this)(*value);
        }
        double lower_bound, upper_bound;
    };
    BetweenVisitor visitor_;
};

} // namespace csr
} // namespace octopus

#endif
