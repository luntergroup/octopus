// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "measure_factory.hpp"

#include <unordered_map>
#include <functional>

#include "measures_fwd.hpp"
#include "exceptions/user_error.hpp"

namespace octopus { namespace csr {

using MeasureMakerMap = std::unordered_map<std::string, std::function<MeasureWrapper()>>;

void init(MeasureMakerMap& measure_makers)
{
    measure_makers[name<AlleleFrequency>()]          = [] () { return make_wrapped_measure<AlleleFrequency>(); };
    measure_makers[name<Depth>()]                    = [] () { return make_wrapped_measure<Depth>(); };
    measure_makers[name<MappingQualityDivergence>()] = [] () { return make_wrapped_measure<MappingQualityDivergence>(); };
    measure_makers[name<MappingQualityZeroCount>()]  = [] () { return make_wrapped_measure<MappingQualityZeroCount>(); };
    measure_makers[name<MeanMappingQuality>()]       = [] () { return make_wrapped_measure<MeanMappingQuality>(); };
    measure_makers[name<ModelPosterior>()]           = [] () { return make_wrapped_measure<ModelPosterior>(); };
    measure_makers[name<Quality>()]                  = [] () { return make_wrapped_measure<Quality>(); };
    measure_makers[name<QualityByDepth>()]           = [] () { return make_wrapped_measure<QualityByDepth>(); };
    measure_makers[name<GenotypeQuality>()]          = [] () { return make_wrapped_measure<GenotypeQuality>(); };
    measure_makers[name<StrandBias>()]               = [] () { return make_wrapped_measure<StrandBias>(); };
    measure_makers[name<GCContent>()]                = [] () { return make_wrapped_measure<GCContent>(); };
    measure_makers[name<FilteredReadFraction>()]     = [] () { return make_wrapped_measure<FilteredReadFraction>(); };
    measure_makers[name<ClippedReadFraction>()]      = [] () { return make_wrapped_measure<ClippedReadFraction>(); };
    measure_makers[name<IsDenovo>()]                 = [] () { return make_wrapped_measure<IsDenovo>(); };
    measure_makers[name<IsSomatic>()]                = [] () { return make_wrapped_measure<IsSomatic>(); };
    measure_makers[name<AmbiguousReadFraction>()]    = [] () { return make_wrapped_measure<AmbiguousReadFraction>(); };
    measure_makers[name<MedianBaseQuality>()]        = [] () { return make_wrapped_measure<MedianBaseQuality>(); };
    measure_makers[name<AltAlleleCount>()]           = [] () { return make_wrapped_measure<AltAlleleCount>(); };
    measure_makers[name<Realignments>()]             = [] () { return make_wrapped_measure<Realignments>(); };
}

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

MeasureWrapper make_measure(const std::string& name)
{
    static MeasureMakerMap measure_makers {};
    if (measure_makers.empty()) {
        init(measure_makers);
    }
    if (measure_makers.count(name) == 0) {
        throw UnknownMeasure {name};
    }
    return measure_makers.at(name)();
}

} // namespace csr
} // namespace octopus
