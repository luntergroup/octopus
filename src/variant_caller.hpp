//
//  variant_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_caller_hpp
#define Octopus_variant_caller_hpp

#include <vector>
#include <string>
#include <iterator>

#include "common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_transform.hpp"
#include "candidate_variant_generator.hpp"

class GenomicRegion;
class Variant;
class VcfRecord;

namespace Octopus
{
    class VariantCaller
    {
    public:
        enum class RefCallType { Positional, Blocked, None };
        
        using ReadMap = Octopus::ReadMap;
        
        VariantCaller() = delete;
        VariantCaller(ReferenceGenome& reference, CandidateVariantGenerator& candidate_generator,
                      RefCallType refcall_type = RefCallType::None);
        virtual ~VariantCaller() = default;
        
        VariantCaller(const VariantCaller&)            = delete;
        VariantCaller& operator=(const VariantCaller&) = delete;
        VariantCaller(VariantCaller&&)                 = delete;
        VariantCaller& operator=(VariantCaller&&)      = delete;
        
        std::string get_details() const;
        size_t num_buffered_reads() const noexcept;
        std::vector<VcfRecord> call_variants(const GenomicRegion& region, ReadMap reads);
        
    protected:
        ReferenceGenome& reference_;
        
        const RefCallType refcall_type_ = RefCallType::Positional;
        
        bool refcalls_requested() const noexcept;
        
    private:
        CandidateVariantGenerator& candidate_generator_;
        
        bool done_calling(const GenomicRegion& region) const noexcept;
        
        virtual std::string do_get_details() const = 0;
        
        //virtual size_t do_num_buffered_reads() const noexcept = 0;
        
        virtual std::vector<VcfRecord> call_variants(const GenomicRegion& region,
                                                     const std::vector<Variant>& candidates,
                                                     const ReadMap& reads) = 0;
    };
    
    std::vector<Allele>
    generate_callable_alleles(const GenomicRegion& region, const std::vector<Variant>& variants,
                              VariantCaller::RefCallType refcall_type, ReferenceGenome& reference);
    
    std::vector<Allele>
    generate_candidate_reference_alleles(const std::vector<Allele>& callable_alleles,
                                         const std::vector<GenomicRegion>& called_regions,
                                         const std::vector<Variant>& candidates,
                                         VariantCaller::RefCallType refcall_type);
    
    template <typename Map>
    void remove_low_posteriors(Map& map, double min_posterior)
    {
        for (auto it = std::begin(map); it != std::end(map);) {
            if (it->second < min_posterior) {
                it = map.erase(it);
            } else {
                ++it;
            }
        }
    }
    
} // namespace Octopus

#endif
