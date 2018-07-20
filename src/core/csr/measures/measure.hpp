// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef measure_hpp
#define measure_hpp

#include <vector>
#include <string>
#include <memory>
#include <utility>
#include <unordered_map>

#include <boost/variant.hpp>
#include <boost/optional.hpp>
#include <boost/any.hpp>

#include "concepts/equitable.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_record.hpp"
#include "../facets/facet.hpp"

namespace octopus { namespace csr {

class Measure : public Equitable<Measure>
{
public:
    using FacetMap = std::unordered_map<std::string, FacetWrapper>;
    using ResultType = boost::variant<double, boost::optional<double>,
                                      std::vector<double>, std::vector<boost::optional<double>>,
                                      int, boost::optional<int>,
                                      std::vector<int>, boost::optional<std::vector<int>>, std::vector<boost::optional<int>>,
                                      std::size_t, boost::optional<std::size_t>,
                                      std::vector<std::size_t>, std::vector<boost::optional<std::size_t>>,
                                      bool,
                                      std::vector<bool>,
                                      boost::any>;
    enum class ResultCardinality { one, num_alleles, num_samples };
    
    Measure() = default;
    
    Measure(const Measure&)            = default;
    Measure& operator=(const Measure&) = default;
    Measure(Measure&&)                 = default;
    Measure& operator=(Measure&&)      = default;
    
    virtual ~Measure() = default;
    
    std::unique_ptr<Measure> clone() const { return do_clone(); }
    
    ResultType evaluate(const VcfRecord& call, const FacetMap& facets) const { return do_evaluate(call, facets); }
    ResultCardinality cardinality() const noexcept { return do_cardinality(); }
    const std::string& name() const { return do_name(); }
    std::string describe() const { return do_describe(); }
    std::vector<std::string> requirements() const { return do_requirements(); }
    std::string serialise(const ResultType& value) const { return do_serialise(value); }
    void annotate(VcfHeader::Builder& header) const;
    void annotate(VcfRecord::Builder& record, const ResultType& value) const;
    friend bool operator==(const Measure& lhs, const Measure& rhs) noexcept
    {
        return lhs.name() == rhs.name() && lhs.is_equal(rhs);
    }
    
private:
    virtual std::unique_ptr<Measure> do_clone() const = 0;
    virtual ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const = 0;
    virtual ResultCardinality do_cardinality() const noexcept = 0;
    virtual const std::string& do_name() const = 0;
    virtual std::string do_describe() const = 0;
    virtual std::vector<std::string> do_requirements() const { return {}; }
    virtual std::string do_serialise(const ResultType& value) const;
    virtual bool is_required_vcf_field() const noexcept { return false; }
    virtual bool is_equal(const Measure& other) const noexcept { return true; }
};

class MeasureWrapper : public Equitable<MeasureWrapper>
{
public:
    MeasureWrapper() = delete;
    
    MeasureWrapper(std::unique_ptr<Measure> measure) : measure_ {std::move(measure)} {}
    
    MeasureWrapper(const MeasureWrapper& other) : measure_ {other.measure_->clone()} {}
    MeasureWrapper& operator=(const MeasureWrapper& other)
    {
        if (&other != this) measure_ = other.measure_->clone();
        return *this;
    }
    MeasureWrapper(MeasureWrapper&&)            = default;
    MeasureWrapper& operator=(MeasureWrapper&&) = default;
    
    ~MeasureWrapper() = default;
    
    const Measure* base() const noexcept { return measure_.get(); }
    auto operator()(const VcfRecord& call) const { return measure_->evaluate(call, {}); }
    auto operator()(const VcfRecord& call, const Measure::FacetMap& facets) const { return measure_->evaluate(call, facets); }
    Measure::ResultCardinality cardinality() const noexcept { return measure_->cardinality(); }
    const std::string& name() const { return measure_->name(); }
    std::string describe() const { return measure_->describe(); }
    std::vector<std::string> requirements() const { return measure_->requirements(); }
    std::string serialise(const Measure::ResultType& value) const { return measure_->serialise(value); }
    void annotate(VcfHeader::Builder& header) const { measure_->annotate(header); }
    void annotate(VcfRecord::Builder& record, const Measure::ResultType& value) const { measure_->annotate(record, value); }
    
private:
    std::unique_ptr<Measure> measure_;
};

inline bool operator==(const MeasureWrapper& lhs, const MeasureWrapper& rhs) noexcept
{
    return *lhs.base() == *rhs.base();
}

template <typename M, typename... Args>
MeasureWrapper make_wrapped_measure(Args&&... args)
{
    return MeasureWrapper {std::make_unique<M>(std::forward<Args>(args)...)};
}

template <typename Measure>
const std::string& name()
{
    return Measure().name();
}

bool is_missing(const Measure::ResultType& value) noexcept;

Measure::ResultType get_sample_value(const Measure::ResultType& value, const MeasureWrapper& measure, std::size_t sample_idx);
std::vector<Measure::ResultType> get_sample_values(const std::vector<Measure::ResultType>& values,
                                                   const std::vector<MeasureWrapper>& measures,
                                                   std::size_t sample_idx);

} // namespace csr
} // namespace octopus

namespace std {

template <> struct hash<octopus::csr::MeasureWrapper>
{
    size_t operator()(const octopus::csr::MeasureWrapper& measure) const
    {
        return hash<string>{}(measure.name());
    }
};

} // namespace std

#endif
