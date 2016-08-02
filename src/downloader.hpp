//
//  downloader.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/06/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__downloader__
#define __Octopus__downloader__

#include <vector>
#include <functional>
#include <memory>

#include "variant_generator.hpp"
#include "variant.hpp"

namespace octopus {

class GenomicRegion;
class ReferenceGenome;
class AlignedRead;

namespace coretools {

class Downloader : public VariantGenerator
{
public:
    Downloader() = delete;
    
    Downloader(const ReferenceGenome& reference, Variant::RegionType::Size max_variant_size = 100);
    
    Downloader(const Downloader&)            = default;
    Downloader& operator=(const Downloader&) = default;
    Downloader(Downloader&&)                 = default;
    Downloader& operator=(Downloader&&)      = default;
    
    ~Downloader() override = default;
    
private:
    std::unique_ptr<VariantGenerator> do_clone() const override;
    
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    
    std::string name() const override;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Variant::RegionType::Size max_variant_size_;
};

} // namespace coretools
} // namespace octopus

#endif /* defined(__Octopus__downloader__) */
