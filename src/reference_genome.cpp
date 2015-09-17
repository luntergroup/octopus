//
//  reference_genome.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "reference_genome.h"

#include <algorithm> // std::find, std::sort, std::remove_copy
#include <iterator>  // std::cbegin, std::cend, std::back_inserter
#include <utility>   // std::move
#include <regex>     // std::regex, std::smatch, std::regex_search
#include <stdexcept>

#include "fasta.h"
#include "caching_fasta.h"
#include "threadsafe_fasta.h"

ReferenceGenome::ReferenceGenome(std::unique_ptr<ReferenceGenomeImpl> impl)
:
impl_ {std::move(impl)},
name_ {impl_->get_reference_name()},
contig_names_(std::move(impl_->get_contig_names())),
contig_sizes_ {}
{
    for (const auto& contig_name : contig_names_) {
        contig_sizes_[contig_name] = impl_->get_contig_size(contig_name);
    }
}

const std::string& ReferenceGenome::get_name() const
{
    return name_;
}

bool ReferenceGenome::has_contig(const std::string& contig_name) const noexcept
{
    return std::find(std::cbegin(contig_names_), std::cend(contig_names_), contig_name) != std::cend(contig_names_);
}

std::size_t ReferenceGenome::num_contigs() const noexcept
{
    return contig_names_.size();
}

const std::vector<std::string>& ReferenceGenome::get_contig_names() const noexcept
{
    return contig_names_;
}

ReferenceGenome::SizeType ReferenceGenome::get_contig_size(const std::string& contig_name) const
{
    if (has_contig(contig_name)) {
        return contig_sizes_.at(contig_name);
    }
    throw std::runtime_error {"contig " + contig_name + " is not in reference genome " + name_};
}

ReferenceGenome::SizeType ReferenceGenome::get_contig_size(const GenomicRegion& region) const
{
    return get_contig_size(region.get_contig_name());
}

GenomicRegion ReferenceGenome::get_contig_region(const std::string& contig_name) const
{
    return GenomicRegion {contig_name, 0, get_contig_size(contig_name)};
}

bool ReferenceGenome::contains_region(const GenomicRegion& region) const noexcept
{
    return region.get_end() <= get_contig_size(region);
}

ReferenceGenome::SequenceType ReferenceGenome::get_sequence(const GenomicRegion& region)
{
    return impl_->get_sequence(region);
}

// non-member functions

ReferenceGenome make_reference(fs::path file_path, std::size_t max_base_pair_cache, bool is_threaded)
{
    if (max_base_pair_cache > 0) {
        return ReferenceGenome {std::make_unique<CachingFasta>(std::move(file_path), max_base_pair_cache)};
    } else if (is_threaded) {
        return ReferenceGenome {std::make_unique<ThreadsafeFasta>(std::move(file_path))};
    } else {
        return ReferenceGenome {std::make_unique<Fasta>(std::move(file_path))};
    }
}

std::vector<GenomicRegion> get_all_contig_regions(const ReferenceGenome& reference)
{
    std::vector<GenomicRegion> result {};
    result.reserve(reference.num_contigs());
    
    for (const auto& contig : reference.get_contig_names()) {
        result.emplace_back(reference.get_contig_region(contig));
    }
    
    std::sort(result.begin(), result.end(), [] (const auto& lhs, const auto& rhs) {
        return size(lhs) < size(rhs);
    });
    
    return result;
}

// Requires reference access to get contig sizes for partially specified regions (e.g. "4")
GenomicRegion parse_region(const std::string& region, const ReferenceGenome& reference)
{
    std::string filtered_region;
    std::remove_copy(region.cbegin(), region.cend(), std::back_inserter(filtered_region), ',');
    
    const static std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};
    std::smatch match;
    
    if (std::regex_search(filtered_region, match, re) && match.size() == 5) {
        auto contig_name = match.str(1);
        GenomicRegion::SizeType begin {}, end {};
        auto the_contig_size = reference.get_contig_size(contig_name);
        
        if (match.str(2).empty()) {
            end = the_contig_size;
        } else {
            begin = static_cast<GenomicRegion::SizeType>(std::stoul(match.str(2)));
            
            if (match.str(3).empty()) {
                end = begin;
            } else if (match.str(4).empty()) {
                end = the_contig_size;
            } else {
                end = static_cast<GenomicRegion::SizeType>(std::stoul(match.str(4)));
            }
            
            if (begin > the_contig_size || end > the_contig_size) {
                throw std::runtime_error {"region " + region + " is larger than contig in " + reference.get_name()};
            }
        }
        
        return GenomicRegion {std::move(contig_name), begin, end};
    }
    
    throw std::runtime_error {"region" + region + " has invalid format"};
}
