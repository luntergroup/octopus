// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef vcf_extractor_hpp
#define vcf_extractor_hpp

#include <vector>
#include <memory>

#include <boost/filesystem.hpp>
#include <boost/optional.hpp>

#include "io/variant/vcf.hpp"
#include "core/types/variant.hpp"
#include "variant_generator.hpp"

namespace octopus {

class GenomicRegion;

namespace coretools {

class VcfExtractor : public VariantGenerator
{
public:
    struct Options
    {
        Variant::MappingDomain::Size max_variant_size = 100;
        bool extract_filtered = false;
        boost::optional<VcfRecord::QualityType> min_quality = boost::none;
    };
    
    VcfExtractor() = delete;
    
    VcfExtractor(std::unique_ptr<const VcfReader> reader);
    VcfExtractor(std::unique_ptr<const VcfReader> reader, Options options);
        
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
    Options options_;
    
    bool is_good(const VcfRecord& record);
};

} // namespace coretools
} // namespace octopus

#endif
