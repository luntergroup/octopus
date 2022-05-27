// Copyright (c) 2015-2021 Daniel Cooke
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
        bool split_complex = false;
    };
    
    VcfExtractor() = delete;
    
    VcfExtractor(std::unique_ptr<VcfReader> reader);
    VcfExtractor(std::unique_ptr<VcfReader> reader, Options options);
    
    VcfExtractor(const VcfExtractor&)            = default;
    VcfExtractor& operator=(const VcfExtractor&) = default;
    VcfExtractor(VcfExtractor&&)                 = default;
    VcfExtractor& operator=(VcfExtractor&&)      = default;
    
    ~VcfExtractor() override = default;
    
private:
    std::unique_ptr<VariantGenerator> do_clone() const override;
    std::vector<Variant> do_generate(const RegionSet& regions, OptionalThreadPool workers) const override;
    std::string name() const override;
    
    mutable std::shared_ptr<VcfReader> reader_;
    Options options_;
    
    std::vector<Variant> fetch_variants(const GenomicRegion& region) const;
    bool is_good(const VcfRecord& record) const;
};

} // namespace coretools
} // namespace octopus

#endif
