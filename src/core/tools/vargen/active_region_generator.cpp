// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "active_region_generator.hpp"

namespace octopus { namespace coretools {

ActiveRegionGenerator::ActiveRegionGenerator(const ReferenceGenome& reference, Options options)
: reference_ {reference}
, options_ {options}
, assembler_name_ {"LocalReassembler"}
, cigar_scanner_name_ {"CigarScanner"}
, assembler_active_region_generator_ {}
{}

void ActiveRegionGenerator::add_generator(const std::string& name)
{
    if (is_assembler(name) && !options_.assemble_all) {
        assembler_active_region_generator_ = AssemblerActiveRegionGenerator {reference_};
    }
}

void ActiveRegionGenerator::add_read(const SampleName& sample, const AlignedRead& read)
{
    if (assembler_active_region_generator_) assembler_active_region_generator_->add(sample, read);
}

std::vector<GenomicRegion> ActiveRegionGenerator::generate(const GenomicRegion& region, const std::string& generator) const
{
    if (using_assembler()) {
        if (is_assembler(generator)) {
            if (assembler_active_region_generator_) {
                return assembler_active_region_generator_->generate(region);
            } else {
                return {region};
            }
        } else if (is_cigar_scanner(generator)) {
            return {region};
        }
    }
    return {region};
}

void ActiveRegionGenerator::clear() noexcept
{
    if (assembler_active_region_generator_) assembler_active_region_generator_->clear();
}

// private methods

bool ActiveRegionGenerator::is_cigar_scanner(const std::string& generator) const noexcept
{
    return generator == cigar_scanner_name_;
}

bool ActiveRegionGenerator::is_assembler(const std::string& generator) const noexcept
{
    return generator == assembler_name_;
}

bool ActiveRegionGenerator::using_assembler() const noexcept
{
    return assembler_active_region_generator_.is_initialized();
}

} // namespace coretools
} // namespace octopus
