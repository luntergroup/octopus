//
//  vcf_extractor.hpp
//  octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__vcf_extractor__
#define __Octopus__vcf_extractor__

#include <vector>
#include <memory>

#include <boost/filesystem.hpp>

#include <io/variant/vcf_reader.hpp>
#include <core/types/variant.hpp>

#include "variant_generator.hpp"

namespace octopus {

class GenomicRegion;

namespace coretools {

class VcfExtractor : public VariantGenerator
{
public:
    VcfExtractor() = delete;
    
    VcfExtractor(std::unique_ptr<const VcfReader> reader);
    
    VcfExtractor(const std::shared_ptr<const VcfReader>& reader);
    
    VcfExtractor(const VcfExtractor&)            = default;
    VcfExtractor& operator=(const VcfExtractor&) = default;
    VcfExtractor(VcfExtractor&&)                 = default;
    VcfExtractor& operator=(VcfExtractor&&)      = default;
    
    ~VcfExtractor() override = default;
    
private:
    std::unique_ptr<VariantGenerator> do_clone() const override;
    
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    
    std::string name() const override;
    
    std::shared_ptr<const VcfReader> reader_;
};

} // namespace coretools
} // namespace octopus

#endif /* defined(__Octopus__vcf_extractor__) */
