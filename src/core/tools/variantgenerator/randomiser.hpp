// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef randomiser_hpp
#define randomiser_hpp

#include <vector>
#include <functional>
#include <memory>

#include "variant_generator.hpp"
#include <basics/aligned_read.hpp>
#include <core/types/variant.hpp>

namespace octopus {
    
class GenomicRegion;
class ReferenceGenome;

namespace coretools {

class Randomiser : public VariantGenerator
{
public:
    Randomiser() = delete;
    
    Randomiser(const ReferenceGenome& reference);
    
    Randomiser(const Randomiser&)            = default;
    Randomiser& operator=(const Randomiser&) = default;
    Randomiser(Randomiser&&)                 = default;
    Randomiser& operator=(Randomiser&&)      = default;
    
    ~Randomiser() override = default;
    
private:
    using VariantGenerator::VectorIterator;
    using VariantGenerator::FlatSetIterator;
    
    std::unique_ptr<VariantGenerator> do_clone() const override;
    
    void do_add_reads(VectorIterator first, VectorIterator last) override;
    void do_add_reads(FlatSetIterator first, FlatSetIterator last) override;
    
    std::vector<Variant> do_generate_variants(const GenomicRegion& region) override;
    
    std::string name() const override;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    AlignedRead::MappingDomain::Size max_read_size_ = 100;
};

} // namespace coretools
} // namespace octopus

#endif /* randomiser_hpp */
