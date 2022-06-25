// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "measure_factory.hpp"

#include <unordered_map>
#include <functional>
#include <algorithm>

#include "measures_fwd.hpp"
#include "exceptions/user_error.hpp"
#include "utils/map_utils.hpp"
#include "utils/string_utils.hpp"
#include "io/variant/vcf_header.hpp"
#include "io/variant/vcf_spec.hpp"

namespace octopus { namespace csr {

using MeasureMakerMap = std::unordered_map<std::string, std::function<MeasureWrapper()>>;

void init(MeasureMakerMap& measure_makers)
{
    measure_makers[name<AlleleDepth>()]                 = [] () { return make_wrapped_measure<AlleleDepth>(); };
    measure_makers[name<AlleleFrequency>()]             = [] () { return make_wrapped_measure<AlleleFrequency>(); };
    measure_makers[name<AlleleFrequencyBias>()]         = [] () { return make_wrapped_measure<AlleleFrequencyBias>(); };
    measure_makers[name<Depth>()]                       = [] () { return make_wrapped_measure<Depth>(); };
    measure_makers[name<MappingQualityDivergence>()]    = [] () { return make_wrapped_measure<MappingQualityDivergence>(); };
    measure_makers[name<MappingQualityZeroCount>()]     = [] () { return make_wrapped_measure<MappingQualityZeroCount>(); };
    measure_makers[name<MeanMappingQuality>()]          = [] () { return make_wrapped_measure<MeanMappingQuality>(); };
    measure_makers[name<ModelPosterior>()]              = [] () { return make_wrapped_measure<ModelPosterior>(); };
    measure_makers[name<ModelPosteriorByDepth>()]       = [] () { return make_wrapped_measure<ModelPosteriorByDepth>(); };
    measure_makers[name<Quality>()]                     = [] () { return make_wrapped_measure<Quality>(); };
    measure_makers[name<QualityByDepth>()]              = [] () { return make_wrapped_measure<QualityByDepth>(); };
    measure_makers[name<GenotypeQuality>()]             = [] () { return make_wrapped_measure<GenotypeQuality>(); };
    measure_makers[name<GenotypeQualityByDepth>()]      = [] () { return make_wrapped_measure<GenotypeQualityByDepth>(); };
    measure_makers[name<StrandBias>()]                  = [] () { return make_wrapped_measure<StrandBias>(); };
    measure_makers[name<GCContent>()]                   = [] () { return make_wrapped_measure<GCContent>(); };
    measure_makers[name<FilteredReadFraction>()]        = [] () { return make_wrapped_measure<FilteredReadFraction>(); };
    measure_makers[name<ClippedReadFraction>()]         = [] () { return make_wrapped_measure<ClippedReadFraction>(); };
    measure_makers[name<IsDenovo>()]                    = [] () { return make_wrapped_measure<IsDenovo>(); };
    measure_makers[name<IsSomatic>()]                   = [] () { return make_wrapped_measure<IsSomatic>(); };
    measure_makers[name<AmbiguousReadFraction>()]       = [] () { return make_wrapped_measure<AmbiguousReadFraction>(); };
    measure_makers[name<MedianBaseQuality>()]           = [] () { return make_wrapped_measure<MedianBaseQuality>(); };
    measure_makers[name<MismatchCount>()]               = [] () { return make_wrapped_measure<MismatchCount>(); };
    measure_makers[name<MismatchFraction>()]            = [] () { return make_wrapped_measure<MismatchFraction>(); };
    measure_makers[name<IsRefcall>()]                   = [] () { return make_wrapped_measure<IsRefcall>(); };
    measure_makers[name<NormalContamination>()]         = [] () { return make_wrapped_measure<NormalContamination>(); };
    measure_makers[name<DeNovoContamination>()]         = [] () { return make_wrapped_measure<DeNovoContamination>(); };
    measure_makers[name<ReadSideBias>()]                = [] () { return make_wrapped_measure<ReadSideBias>(); };
    measure_makers[name<AltAlleleCount>()]              = [] () { return make_wrapped_measure<AltAlleleCount>(); };
    measure_makers[name<STRLength>()]                   = [] () { return make_wrapped_measure<STRLength>(); };
    measure_makers[name<STRPeriod>()]                   = [] () { return make_wrapped_measure<STRPeriod>(); };
    measure_makers[name<PosteriorProbability>()]        = [] () { return make_wrapped_measure<PosteriorProbability>(); };
    measure_makers[name<PosteriorProbabilityByDepth>()] = [] () { return make_wrapped_measure<PosteriorProbabilityByDepth>(); };
    measure_makers[name<ClassificationConfidence>()]    = [] () { return make_wrapped_measure<ClassificationConfidence>(); };
    measure_makers[name<SomaticHaplotypeCount>()]       = [] () { return make_wrapped_measure<SomaticHaplotypeCount>(); };
    measure_makers[name<MedianSomaticMappingQuality>()] = [] () { return make_wrapped_measure<MedianSomaticMappingQuality>(); };
    measure_makers[name<StrandDisequilibrium>()]        = [] () { return make_wrapped_measure<StrandDisequilibrium>(); };
    measure_makers[name<SupplementaryFraction>()]       = [] () { return make_wrapped_measure<SupplementaryFraction>(); };
    measure_makers[name<MisalignedReadCount>()]         = [] () { return make_wrapped_measure<MisalignedReadCount>(); };
    measure_makers[name<ReadTailBias>()]                = [] () { return make_wrapped_measure<ReadTailBias>(); };
    measure_makers[name<ReadEndBias>()]                 = [] () { return make_wrapped_measure<ReadEndBias>(); };
    measure_makers[name<VariantLength>()]               = [] () { return make_wrapped_measure<VariantLength>(); };
    measure_makers[name<BaseMismatchCount>()]           = [] () { return make_wrapped_measure<BaseMismatchCount>(); };
    measure_makers[name<BaseMismatchFraction>()]        = [] () { return make_wrapped_measure<BaseMismatchFraction>(); };
    measure_makers[name<BaseMismatchQuality>()]         = [] () { return make_wrapped_measure<BaseMismatchQuality>(); };
    measure_makers[name<AssignedDepth>()]               = [] () { return make_wrapped_measure<AssignedDepth>(); };
    measure_makers[name<DuplicateConcordance>()]        = [] () { return make_wrapped_measure<DuplicateConcordance>(); };
    measure_makers[name<DuplicateAlleleDepth>()]        = [] () { return make_wrapped_measure<DuplicateAlleleDepth>(); };
    measure_makers[name<DuplicateAlleleFraction>()]     = [] () { return make_wrapped_measure<DuplicateAlleleFraction>(); };
    measure_makers[name<ErrorRate>()]                   = [] () { return make_wrapped_measure<ErrorRate>(); };
    measure_makers[name<ErrorRateStdev>()]              = [] () { return make_wrapped_measure<ErrorRateStdev>(); };
    measure_makers[name<IsTransversion>()]              = [] () { return make_wrapped_measure<IsTransversion>(); };
    measure_makers[name<PhaseLength>()]                 = [] () { return make_wrapped_measure<PhaseLength>(); };
    measure_makers[name<MaxReadLength>()]               = [] () { return make_wrapped_measure<MaxReadLength>(); };
    measure_makers[name<AlleleMappingQuality>()]        = [] () { return make_wrapped_measure<AlleleMappingQuality>(); };
    measure_makers[name<MeanLikelihood>()]              = [] () { return make_wrapped_measure<MeanLikelihood>(); };
    measure_makers[name<PhylogenyPosterior>()]          = [] () { return make_wrapped_measure<PhylogenyPosterior>(); };
}

class BadParameterList : public UserError
{
    std::string name_;
    std::string do_where() const override { return "make_measure"; }
    std::string do_why() const override
    {
        return name_ + std::string {" does not list parameters in required format"};
    }
    std::string do_help() const override
    {
        return "Format is MEASURE[param1,param2,...,paramn]";
    }
public:
    BadParameterList(std::string name) : name_ {std::move(name)} {}
};

class UnknownMeasure : public UserError
{
    std::string name_;
    std::string do_where() const override { return "make_measure"; }
    std::string do_why() const override
    {
        return name_ + std::string {" is not a valid measure"};
    }
    std::string do_help() const override
    {
        return "See the documentation for valid measures";
    }
public:
    UnknownMeasure(std::string name) : name_ {std::move(name)} {}
};

bool has_parameters(const std::string& name) noexcept
{
    assert(!name.empty());
    return name.back() == ']';
}

std::vector<std::string> parse_parameters(std::string& name)
{
    assert(has_parameters(name));
    const auto parameter_list_begin = name.find('[');
    if (parameter_list_begin == std::string::npos) {
        throw BadParameterList {name};
    }
    auto parameter_list_str = name.substr(parameter_list_begin + 1, name.size() - parameter_list_begin - 2);
    parameter_list_str.erase(std::remove(parameter_list_str.begin(), parameter_list_str.end(), ' '), parameter_list_str.end());
    auto result = utils::split(parameter_list_str, ',');
    name.erase(parameter_list_begin);
    return result;
}

MeasureWrapper make_measure(std::string name)
{
    static MeasureMakerMap measure_makers {};
    if (measure_makers.empty()) {
        init(measure_makers);
    }
    std::vector<std::string> params {};
    if (has_parameters(name)) {
        params = parse_parameters(name);
    }
    if (measure_makers.count(name) == 0) {
        throw UnknownMeasure {name};
    }
    auto result = measure_makers.at(name)();
    if (!params.empty()) {
        result.set_parameters(std::move(params));
    }
    return result;
}

std::vector<MeasureWrapper> make_measures(std::vector<std::string> names)
{
    std::vector<MeasureWrapper> result {};
    result.reserve(names.size());
    std::transform(std::cbegin(names), std::cend(names), std::back_inserter(result), make_measure);
    return result;
}

std::vector<std::string> get_all_measure_names()
{
    static MeasureMakerMap measure_makers {};
    if (measure_makers.empty()) {
        init(measure_makers);
    }
    return extract_sorted_keys(measure_makers);
}

VcfHeader make_header(const std::vector<MeasureWrapper>& measures)
{
    VcfHeader::Builder builder {};
    for (const auto& measure : measures) {
        measure.annotate(builder);
    }
    return builder.build_once();
}

void print_help(const std::vector<MeasureWrapper>& measures, std::ostream& os)
{
    const auto header = make_header(measures); 
    os << "Name\tKind\tNumber\tType\tDescription" << std::endl;
    for (const auto& measure : measures) {
        os << measure.name() << '\t';
        using namespace vcfspec::header;
        const auto kind = header.has(format, measure.name()) ? format : info;
        os << kind << '\t';
        os << get_id_field_value(header, kind, measure.name(), meta::struc::number) << '\t';
        os << get_id_field_value(header, kind, measure.name(), meta::struc::type) << '\t';
        os << get_id_field_value(header, kind, measure.name(), meta::struc::description) << '\t';
        os << std::endl;
    }
}

void print_all_measures_help(std::ostream& os)
{
    const auto measures = make_measures(get_all_measure_names());
    print_help(measures, os);
}

} // namespace csr
} // namespace octopus
