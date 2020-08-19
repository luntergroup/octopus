// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_generator_builder_hpp
#define variant_generator_builder_hpp

#include <deque>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "cigar_scanner.hpp"
#include "local_reassembler.hpp"
#include "vcf_extractor.hpp"
#include "repeat_scanner.hpp"
#include "downloader.hpp"
#include "randomiser.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/variant/vcf_reader.hpp"
#include "active_region_generator.hpp"

namespace octopus {
namespace coretools {

class VariantGenerator;

class VariantGeneratorBuilder
{
public:
    VariantGeneratorBuilder() = default;
    
    VariantGeneratorBuilder(const VariantGeneratorBuilder&)            = default;
    VariantGeneratorBuilder& operator=(const VariantGeneratorBuilder&) = default;
    VariantGeneratorBuilder(VariantGeneratorBuilder&&)                 = default;
    VariantGeneratorBuilder& operator=(VariantGeneratorBuilder&&)      = default;
    
    ~VariantGeneratorBuilder() = default;
    
    VariantGeneratorBuilder& set_cigar_scanner(CigarScanner::Options options);
    VariantGeneratorBuilder& set_local_reassembler(LocalReassembler::Options options);
    VariantGeneratorBuilder& add_vcf_extractor(boost::filesystem::path reader,
                                               VcfExtractor::Options options = VcfExtractor::Options {});
    VariantGeneratorBuilder& set_repeat_scanner(RepeatScanner::Options options);
    VariantGeneratorBuilder& add_downloader(Downloader::Options options = Downloader::Options {});
    VariantGeneratorBuilder& add_randomiser(Randomiser::Options options = Randomiser::Options {});
    VariantGeneratorBuilder& set_active_region_generator(ActiveRegionGenerator::Options options = ActiveRegionGenerator::Options {});
    
    VariantGenerator build(const ReferenceGenome& reference) const;

private:
    struct VcfExtractorPacket
    {
        boost::filesystem::path file;
        VcfExtractor::Options options;
    };
    
    boost::optional<CigarScanner::Options> cigar_scanner_;
    boost::optional<LocalReassembler::Options> local_reassembler_;
    std::deque<VcfExtractorPacket> vcf_extractors_;
    boost::optional<RepeatScanner::Options> repeat_scanner_;
    std::deque<Downloader::Options> downloaders_;
    std::deque<Randomiser::Options> randomisers_;
    ActiveRegionGenerator::Options active_region_generator_;
};

} // namespace coretools

using coretools::VariantGeneratorBuilder;

} // namespace octopus

#endif
