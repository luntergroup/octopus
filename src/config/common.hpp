// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef Octopus_common_hpp
#define Octopus_common_hpp

#include <string>
#include <vector>
#include <cstdint>

#include <boost/optional.hpp>

#include <basics/genomic_region.hpp>
#include <basics/aligned_read.hpp>

#include <io/read/read_manager.hpp>
#include <containers/mappable_flat_set.hpp>
#include <containers/mappable_flat_multi_set.hpp>
#include <containers/mappable_map.hpp>
#include <logging/logging.hpp>

namespace octopus {

extern bool DEBUG_MODE;
extern bool TRACE_MODE;

namespace info
{
    const static std::string VERSION {"1.0"};
    const static std::string BUG_EMAIL {"dcooke@well.ox.ac.uk"};
    const static std::vector<std::string> AUTHORS {"Daniel Cooke"};
    const static std::string COPYRIGHT_NOTICE {"Copyright (c) 2016 University of Oxford"};
}

using SampleName = std::string;

using ContigName = GenomicRegion::ContigName;

using InputRegionMap = MappableSetMap<ContigName, GenomicRegion>;

using ReadContainer = MappableFlatMultiSet<AlignedRead>;
using ReadMap       = MappableMap<SampleName, AlignedRead>;

void log_program_startup(); // Always uses InfoLogger

template <typename Log>
void log_program_end(Log& log)
{
    log << "------------------------------------------------------------------------";
}

inline void log_program_end()
{
    logging::InfoLogger log {};
    log_program_end(log);
}

namespace logging {
    boost::optional<DebugLogger> get_debug_log();
    boost::optional<TraceLogger> get_trace_log();
}

} // namespace octopus

#endif
