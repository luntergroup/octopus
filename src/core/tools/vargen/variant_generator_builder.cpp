// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variant_generator_builder.hpp"

#include <memory>
#include <utility>

#include "variant_generator.hpp"

namespace octopus { namespace coretools {

VariantGeneratorBuilder&
VariantGeneratorBuilder::set_cigar_scanner(CigarScanner::Options options)
{
    cigar_scanner_ = std::move(options);
    return *this;
}

VariantGeneratorBuilder&
VariantGeneratorBuilder::set_dynamic_cigar_scanner(DynamicCigarScanner::Options options)
{
    dynamic_cigar_scanner_ = std::move(options);
    return *this;
}

VariantGeneratorBuilder&
VariantGeneratorBuilder::set_local_reassembler(LocalReassembler::Options options)
{
    local_reassembler_ = std::move(options);
    return *this;
}

VariantGeneratorBuilder&
VariantGeneratorBuilder::add_vcf_extractor(boost::filesystem::path file, VcfExtractor::Options options)
{
    vcf_extractors_.push_back({std::move(file), std::move(options)});
    return *this;
}

VariantGeneratorBuilder&
VariantGeneratorBuilder::add_downloader(Downloader::Options options)
{
    downloaders_.push_back(std::move(options));
    return *this;
}

VariantGeneratorBuilder&
VariantGeneratorBuilder::add_randomiser(Randomiser::Options options)
{
    randomisers_.push_back(std::move(options));
    return *this;
}

VariantGenerator VariantGeneratorBuilder::build(const ReferenceGenome& reference) const
{
    VariantGenerator result {};
    
    if (cigar_scanner_) {
        result.add(std::make_unique<CigarScanner>(reference, *cigar_scanner_));
    }
    if (dynamic_cigar_scanner_) {
        result.add(std::make_unique<DynamicCigarScanner>(reference, *dynamic_cigar_scanner_));
    }
    if (local_reassembler_) {
        result.add(std::make_unique<LocalReassembler>(reference, *local_reassembler_));
    }
    for (auto packet : vcf_extractors_) {
        result.add(std::make_unique<VcfExtractor>(std::make_unique<VcfReader>(packet.file), packet.options));
    }
    for (auto options : downloaders_) {
        result.add(std::make_unique<Downloader>(reference, options));
    }
    for (auto options : randomisers_) {
        result.add(std::make_unique<Randomiser>(reference, options));
    }
    
    return result;
}
    
} // namespace coretools
} // namespace octopus
