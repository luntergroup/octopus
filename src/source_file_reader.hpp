//
//  source_file_reader.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__source_file_reader__
#define __Octopus__source_file_reader__

#include <vector>
#include <memory>

#include <boost/filesystem.hpp>

#include "candidate_variant_generator.hpp"
#include "vcf_reader.hpp"
#include "variant.hpp"

class GenomicRegion;

namespace octopus { namespace core { namespace generators
{
class SourceFileReader : public CandidateVariantGenerator
{
public:
    SourceFileReader() = delete;
    
    SourceFileReader(boost::filesystem::path path);
    SourceFileReader(std::unique_ptr<const VcfReader> reader);
    SourceFileReader(const std::shared_ptr<const VcfReader>& reader);
    
    SourceFileReader(const SourceFileReader&)            = default;
    SourceFileReader& operator=(const SourceFileReader&) = default;
    SourceFileReader(SourceFileReader&&)                 = default;
    SourceFileReader& operator=(SourceFileReader&&)      = default;
    
    ~SourceFileReader() override = default;
    
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
    
private:
    std::shared_ptr<const VcfReader> reader_;
};
} // namespace generators
} // namespace core
} // namespace octopus

#endif /* defined(__Octopus__source_file_reader__) */
