// Copyright (c) 2015-2020 Daniel Cooke
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

#include "concepts/equitable.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_record.hpp"
#include "exceptions/user_error.hpp"
#include "../facets/facet.hpp"

namespace octopus { namespace csr {

class Measure : public Equitable<Measure>
{
public:
    using FacetMap = std::unordered_map<std::string, FacetWrapper>;

    using ValueType = boost::variant<bool, int, std::size_t, double>;
    template <typename T>
    using Array = std::vector<T>;
    template <typename T>
    using Optional = boost::optional<T>;
    using ResultType = boost::variant<
            // one value per call
            ValueType,
            Optional<ValueType>,
            // one value per (alt) allele or sample
            Array<ValueType>, 
            Array<Optional<ValueType>>,
            Optional<Array<ValueType>>,
            Optional<Array<Optional<ValueType>>>,
            // one value per sample and (alt) allele
            Array<Array<ValueType>>,
            Array<Array<Optional<ValueType>>>,
            Array<Optional<Array<ValueType>>>,
            Array<Optional<Array<Optional<ValueType>>>>
        >;

    enum class ResultCardinality { one, alleles, alt_alleles, samples, samples_and_alleles, samples_and_alt_alleles };

    enum class Aggregator { sum, min, min_tail, max, max_tail, mean };
    
    Measure() = default;
    
    Measure(const Measure&)            = default;
    Measure& operator=(const Measure&) = default;
    Measure(Measure&&)                 = default;
    Measure& operator=(Measure&&)      = default;
    
    virtual ~Measure() = default;
    
    std::unique_ptr<Measure> clone() const { return do_clone(); }
    
    void set_parameters(std::vector<std::string> params) { do_set_parameters(std::move(params)); }
    std::vector<std::string> parameters() const { return do_parameters(); }
    ResultType evaluate(const VcfRecord& call, const FacetMap& facets) const { return do_evaluate(call, facets); }
    ResultCardinality cardinality() const noexcept { return do_cardinality(); }
    const std::string& name() const { return do_name(); }
    std::string describe() const { return do_describe(); }
    std::vector<std::string> requirements() const { return do_requirements(); }
    std::string serialise(const ResultType& value) const { return do_serialise(value); }
    void annotate(VcfHeader::Builder& header, bool aggregate_alleles) const;
    void annotate(VcfRecord::Builder& record, const ResultType& value, const VcfHeader& header, bool aggregate_alleles) const;
    boost::optional<Aggregator> aggregator() const noexcept { return do_aggregator(); }
    friend bool operator==(const Measure& lhs, const Measure& rhs) noexcept
    {
        return lhs.name() == rhs.name() && lhs.is_equal(rhs);
    }
    
private:
    virtual std::unique_ptr<Measure> do_clone() const = 0;
    virtual void do_set_parameters(std::vector<std::string> params);
    virtual std::vector<std::string> do_parameters() const { return {}; }
    virtual ValueType get_value_type() const = 0;
    virtual ResultType do_evaluate(const VcfRecord& call, const FacetMap& facets) const = 0;
    virtual ResultCardinality do_cardinality() const noexcept = 0;
    virtual const std::string& do_name() const = 0;
    virtual std::string do_describe() const = 0;
    virtual std::vector<std::string> do_requirements() const { return {}; }
    virtual std::string do_serialise(const ResultType& value) const;
    virtual bool is_required_vcf_field() const noexcept { return false; }
    virtual bool is_equal(const Measure& other) const noexcept { return true; }
    virtual boost::optional<Aggregator> do_aggregator() const noexcept { return boost::none; }
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
    
    Measure* base() noexcept { return measure_.get(); }
    const Measure* base() const noexcept { return measure_.get(); }
    void set_parameters(std::vector<std::string> params) { measure_->set_parameters(std::move(params)); }
    std::vector<std::string> parameters() const { return measure_->parameters(); }
    auto operator()(const VcfRecord& call) const { return measure_->evaluate(call, {}); }
    auto operator()(const VcfRecord& call, const Measure::FacetMap& facets) const { return measure_->evaluate(call, facets); }
    Measure::ResultCardinality cardinality() const noexcept { return measure_->cardinality(); }
    const std::string& name() const { return measure_->name(); }
    std::string describe() const { return measure_->describe(); }
    std::vector<std::string> requirements() const { return measure_->requirements(); }
    std::string serialise(const Measure::ResultType& value) const { return measure_->serialise(value); }
    void annotate(VcfHeader::Builder& header, bool aggregate_alleles = false) const { measure_->annotate(header, aggregate_alleles); }
    void annotate(VcfRecord::Builder& record, const Measure::ResultType& value, const VcfHeader& header, bool aggregate_alleles = false) const {
         measure_->annotate(record, value, header, aggregate_alleles); }
    boost::optional<Measure::Aggregator> aggregator() const noexcept { return measure_->aggregator(); }
    
private:
    std::unique_ptr<Measure> measure_;
};

class BadMeasureParameters : public UserError
{
    std::string name_, reason_;
    std::string do_where() const override { return "Measure::do_set_parameters"; }
    std::string do_why() const override
    {
        return name_ + " " + reason_;
    }
    std::string do_help() const override
    {
        return "Check the documentation for measures and their parameters";
    }
public:
    BadMeasureParameters(std::string name)
    : name_ {std::move(name)}
    , reason_ {"does not have any parameters"}
    {}
    BadMeasureParameters(std::string name, std::string reason)
    : name_ {std::move(name)}
    , reason_ {std::move(reason)}
    {
        if (!reason_.empty()) {
            auto first_non_space_pos = reason_.find_first_not_of(" ");
            reason_.erase(0, first_non_space_pos);
        }
    }
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

std::string long_name(const Measure& measure);
std::string long_name(const MeasureWrapper& measure);

bool is_missing(const Measure::ResultType& value) noexcept;

template <typename T>
T get_value_type(const Measure::ResultType& value)
{
    return boost::get<T>(boost::get<Measure::ValueType>(value));
}

std::vector<std::string> get_all_requirements(const std::vector<MeasureWrapper>& measures);

Measure::ResultType
 get_sample_value(const Measure::ResultType& value,
                  const MeasureWrapper& measure,
                  std::size_t sample_idx,
                  bool aggregate = false);
std::vector<Measure::ResultType>
 get_sample_values(const std::vector<Measure::ResultType>& values,
                   const std::vector<MeasureWrapper>& measures, 
                   std::size_t sample_idx,
                   bool aggregate = false);

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
