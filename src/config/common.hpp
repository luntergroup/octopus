// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef common_hpp
#define common_hpp

#include <string>
#include <vector>
#include <cstdint>

#include <boost/optional.hpp>

#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "basics/aligned_template.hpp"
#include "io/read/read_manager.hpp"
#include "containers/mappable_flat_set.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "containers/mappable_map.hpp"
#include "logging/logging.hpp"

namespace octopus {

extern bool DEBUG_MODE;
extern bool TRACE_MODE;

using SampleName = std::string;

using ContigName = GenomicRegion::ContigName;

using InputRegionMap = MappableSetMap<ContigName, GenomicRegion>;

using ReadContainer = MappableFlatMultiSet<AlignedRead>;
using ReadMap       = MappableMap<SampleName, AlignedRead>;

using TemplateContainer = MappableFlatMultiSet<AlignedTemplate>;
using TemplateMap       = MappableMap<SampleName, AlignedTemplate>;

enum class ExecutionPolicy { seq, par, par_vec }; // To match Parallelism TS

namespace logging {
    boost::optional<DebugLogger> get_debug_log();
    boost::optional<TraceLogger> get_trace_log();
}

} // namespace octopus

#endif
