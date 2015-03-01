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
#include <stdexcept>
#include <memory> // std::unique_ptr
#include <algorithm> // std::find
#include <iterator> // std::cbegin etc
#include <regex>

#include "genomic_region.h"
#include "reference_genome_impl.h"

class ReferenceGenome
{
public:
    using SequenceType = std::string;
    using SizeType     = IReferenceGenomeImpl::SizeType;
    
    ReferenceGenome() = delete;
    explicit ReferenceGenome(std::unique_ptr<IReferenceGenomeImpl> the_reference_implementation);
    
    ReferenceGenome(const ReferenceGenome&)            = delete;
    ReferenceGenome& operator=(const ReferenceGenome&) = delete;
    ReferenceGenome(ReferenceGenome&&)                 = default;
    ReferenceGenome& operator=(ReferenceGenome&&)      = default;
    
    const std::string& get_name() const;
    bool has_contig(const std::string& contig_name) const noexcept;
    const std::vector<std::string>& get_contig_names() const noexcept;
    SizeType get_contig_size(const std::string& contig_name) const;
    SizeType get_contig_size(const GenomicRegion& a_region) const;
    GenomicRegion get_contig_region(const std::string& contig_name) const;
    bool contains_region(const GenomicRegion& a_region) const noexcept;
    SequenceType get_sequence(const GenomicRegion& a_region);
    
private:
    std::unique_ptr<IReferenceGenomeImpl> the_reference_implementation_;
    std::string name_;
    std::vector<std::string> contig_names_;
    std::unordered_map<std::string, SizeType> contig_sizes_;
};

inline
ReferenceGenome::ReferenceGenome(std::unique_ptr<IReferenceGenomeImpl> the_reference_implementation)
:the_reference_implementation_ {std::move(the_reference_implementation)},
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

inline ReferenceGenome::SizeType ReferenceGenome::get_contig_size(const std::string& contig_name) const
{
    if (has_contig(contig_name)) {
        return contig_sizes_.at(contig_name);
    }
    throw std::runtime_error {"Contig is not in reference genome"};
}

inline ReferenceGenome::SizeType ReferenceGenome::get_contig_size(const GenomicRegion& a_region) const
{
    return get_contig_size(a_region.get_contig_name());
}

inline GenomicRegion ReferenceGenome::get_contig_region(const std::string& contig_name) const
{
    return GenomicRegion(contig_name, 0, get_contig_size(contig_name));
}

inline bool ReferenceGenome::contains_region(const GenomicRegion& a_region) const noexcept
{
    return a_region.get_end() <= get_contig_size(a_region);
}

inline ReferenceGenome::SequenceType ReferenceGenome::get_sequence(const GenomicRegion& a_region)
{
    return the_reference_implementation_->get_sequence(a_region);
}

inline GenomicRegion parse_region(const std::string& a_region, const ReferenceGenome& the_reference)
{
    const static std::regex re {"([^:]+)(?::(\\d+)-?(\\d*))?"};
    std::smatch match;
    if (std::regex_search(a_region, match, re) && match.size() == 4) {
        auto contig_name = match.str(1);
        GenomicRegion::SizeType begin {}, end {};
        auto the_contig_size = the_reference.get_contig_size(contig_name);
        if (match.str(2).empty()) {
            end = the_contig_size;
        } else {
            begin = static_cast<GenomicRegion::SizeType>(std::stoul(match.str(2)));
            
            if (match.str(3).empty()) {
                end = begin;
            } else {
                end = static_cast<GenomicRegion::SizeType>(std::stoul(match.str(3)));
            }
            
            if (begin > the_contig_size || end > the_contig_size) {
                throw std::runtime_error {"Region out of bounds"};
            }
        }
        return GenomicRegion {std::move(contig_name), begin, end};
    }
    throw std::runtime_error {"Invalid region format"};
}

#endif /* defined(__Octopus__reference_genome__) */
