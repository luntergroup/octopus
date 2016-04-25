//
//  vcf_record_factory.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef vcf_record_factory_hpp
#define vcf_record_factory_hpp

#include <vector>
#include <memory>

#include "common.hpp"
#include "reference_genome.hpp"
#include "vcf_record.hpp"
#include "call.hpp"

namespace Octopus
{
    class VcfRecordFactory
    {
    public:
        VcfRecordFactory() = delete;
        
        explicit VcfRecordFactory(const ReferenceGenome& reference, const ReadMap& reads,
                                  std::vector<SampleIdType> samples, bool sites_only);
        
        ~VcfRecordFactory() = default;
        
        VcfRecordFactory(const VcfRecordFactory&)            = default;
        VcfRecordFactory& operator=(const VcfRecordFactory&) = default;
        VcfRecordFactory(VcfRecordFactory&&)                 = default;
        VcfRecordFactory& operator=(VcfRecordFactory&&)      = default;
        
        std::vector<VcfRecord> make(std::vector<std::unique_ptr<Call>>&& calls) const;
        
    private:
        const ReferenceGenome& reference_;
        const ReadMap& reads_;
        
        std::vector<SampleIdType> samples_;
        
        bool sites_only_;
        
        VcfRecord make(std::unique_ptr<Call> call) const;
        //VcfRecord make_single(const std::vector<std::unique_ptr<Call>>& calls) const;
    };
} // namespace Octopus

#endif /* vcf_record_factory_hpp */
