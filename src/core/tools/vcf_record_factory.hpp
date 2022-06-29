// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef vcf_record_factory_hpp
#define vcf_record_factory_hpp

#include <vector>
#include <memory>

#include "config/common.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/variant/vcf_record.hpp"
#include "core/types/calls/call.hpp"
#include "core/types/calls/call_wrapper.hpp"

namespace octopus {

class VcfRecordFactory
{
public:
    VcfRecordFactory() = delete;
    
    VcfRecordFactory(const ReferenceGenome& reference, const ReadMap& reads,
                     std::vector<SampleName> samples, bool sites_only);
    
    VcfRecordFactory(const VcfRecordFactory&)            = default;
    VcfRecordFactory& operator=(const VcfRecordFactory&) = delete;
    VcfRecordFactory(VcfRecordFactory&&)                 = default;
    VcfRecordFactory& operator=(VcfRecordFactory&&)      = delete;
    
    ~VcfRecordFactory() = default;
    
    std::vector<VcfRecord> make(std::vector<CallWrapper>&& calls) const;
    
private:
    const ReferenceGenome& reference_;
    const ReadMap& reads_;
    std::vector<SampleName> samples_;
    bool sites_only_;
    double max_qual = 100'000;
    bool skip_inconsistent_ploidy = true;
    
    VcfRecord make(std::unique_ptr<Call> call) const;
    VcfRecord make_segment(std::vector<std::unique_ptr<Call>>&& calls) const;
};

} // namespace octopus

#endif
