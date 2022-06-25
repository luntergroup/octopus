// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef downloader_hpp
#define downloader_hpp

#include <vector>
#include <functional>
#include <memory>

#include "core/types/variant.hpp"
#include "variant_generator.hpp"

namespace octopus {

class GenomicRegion;
class ReferenceGenome;
class AlignedRead;

namespace coretools {

class Downloader : public VariantGenerator
{
public:
    struct Options
    {
        Variant::MappingDomain::Size max_variant_size = 100;
    };
    
    Downloader() = delete;
    
    Downloader(const ReferenceGenome& reference, Options options);
    
    Downloader(const Downloader&)            = default;
    Downloader& operator=(const Downloader&) = default;
    Downloader(Downloader&&)                 = default;
    Downloader& operator=(Downloader&&)      = default;
    
    ~Downloader() override = default;
    
private:
    std::reference_wrapper<const ReferenceGenome> reference_;
    Options options_;
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    std::vector<Variant> do_generate(const RegionSet& regions, OptionalThreadPool workers) const override;
    std::string name() const override;
};

} // namespace coretools
} // namespace octopus

#endif
