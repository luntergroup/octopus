//
//  reference_genome.h
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__reference_genome__
#define __Octopus__reference_genome__

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <algorithm>
#include <exception>
#include <memory>
#include <iterator>

#include "genomic_region.h"
#include "reference_genome_implementor.h"

using std::uint_fast32_t;

class ReferenceGenome
{
public:
    ReferenceGenome() = delete;
    ReferenceGenome(std::unique_ptr<IReferenceGenomeImplementor> the_reference_implementation);
    
    ReferenceGenome(const ReferenceGenome&)            = delete;
    ReferenceGenome& operator=(const ReferenceGenome&) = delete;
    ReferenceGenome(ReferenceGenome&&)                 = default;
    ReferenceGenome& operator=(ReferenceGenome&&)      = default;
    
    const std::string& get_name() const;
    bool has_contig(const std::string& contig_name) const noexcept;
    const std::vector<std::string>& get_contig_names() const noexcept;
    uint_fast32_t get_contig_size(const std::string& contig_name) const;
    uint_fast32_t get_contig_size(const GenomicRegion& a_region) const;
    GenomicRegion get_contig_region(const std::string& contig_name) const;
    bool contains_region(const GenomicRegion& a_region) const;
    std::string get_sequence(const GenomicRegion& a_region);
    
private:
    std::unique_ptr<IReferenceGenomeImplementor> the_reference_implementation_;
    std::string name_;
    std::vector<std::string> contig_names_;
    std::unordered_map<std::string, uint_fast32_t> contig_sizes_;
};

class UnknownContig : public std::exception
{
    virtual const char* what() const throw()
    {
        return "Contig is not in reference genome";
    }
};

inline
ReferenceGenome::ReferenceGenome(std::unique_ptr<IReferenceGenomeImplementor> the_reference_implementation)
:   the_reference_implementation_ {std::move(the_reference_implementation)},
    name_ {the_reference_implementation_->get_reference_name()},
    contig_names_(std::move(the_reference_implementation_->get_contig_names())),
    contig_sizes_ {}
{
    for (const auto& contig_name : contig_names_) {
        contig_sizes_[contig_name] = the_reference_implementation_->get_contig_size(contig_name);
    }
}

inline const std::string& ReferenceGenome::get_name() const
{
    return name_;
}

inline bool ReferenceGenome::has_contig(const std::string& contig_name) const noexcept
{
    return std::find(std::cbegin(contig_names_), std::cend(contig_names_), contig_name) != std::cend(contig_names_);
}

inline const std::vector<std::string>& ReferenceGenome::get_contig_names() const noexcept
{
    return contig_names_;
}

inline uint_fast32_t ReferenceGenome::get_contig_size(const std::string& contig_name) const
{
    if (has_contig(contig_name)) {
        return contig_sizes_.at(contig_name);
    }
    throw UnknownContig {};
}

inline uint_fast32_t ReferenceGenome::get_contig_size(const GenomicRegion& a_region) const
{
    return get_contig_size(a_region.get_contig_name());
}

inline GenomicRegion ReferenceGenome::get_contig_region(const std::string& contig_name) const
{
    return GenomicRegion(contig_name, 0, get_contig_size(contig_name));
}

inline bool ReferenceGenome::contains_region(const GenomicRegion& a_region) const
{
    return a_region.get_end() <= get_contig_size(a_region);
}

inline std::string ReferenceGenome::get_sequence(const GenomicRegion& a_region)
{
    return the_reference_implementation_->get_sequence(a_region);
}

#endif /* defined(__Octopus__reference_genome__) */
