// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variant_generator_builder_hpp
#define variant_generator_builder_hpp

#include <deque>

#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "cigar_scanner.hpp"
#include "dynamic_cigar_scanner.hpp"
#include "local_reassembler.hpp"
#include "vcf_extractor.hpp"
#include "downloader.hpp"
#include "randomiser.hpp"
#include "io/reference/reference_genome.hpp"
#include "io/variant/vcf_reader.hpp"

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
    VariantGeneratorBuilder& set_dynamic_cigar_scanner(DynamicCigarScanner::Options options);
    VariantGeneratorBuilder& set_local_reassembler(LocalReassembler::Options options);
    VariantGeneratorBuilder& add_vcf_extractor(boost::filesystem::path reader,
                                               VcfExtractor::Options options = VcfExtractor::Options {});
    VariantGeneratorBuilder& add_downloader(Downloader::Options options = Downloader::Options {});
    VariantGeneratorBuilder& add_randomiser(Randomiser::Options options = Randomiser::Options {});
    
    VariantGenerator build(const ReferenceGenome& reference) const;

private:
    struct VcfExtractorPacket
    {
        boost::filesystem::path file;
        VcfExtractor::Options options;
    };
    
    boost::optional<CigarScanner::Options> cigar_scanner_;
    boost::optional<DynamicCigarScanner::Options> dynamic_cigar_scanner_;
    boost::optional<LocalReassembler::Options> local_reassembler_;
    std::deque<VcfExtractorPacket> vcf_extractors_;
    std::deque<Downloader::Options> downloaders_;
    std::deque<Randomiser::Options> randomisers_;
};

} // namespace coretools

using coretools::VariantGeneratorBuilder;

} // namespace octopus

#endif
