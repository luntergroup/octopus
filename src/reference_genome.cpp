//
//  reference_genome.cpp
//  Octopus
//
//  Created by Daniel Cooke on 18/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "reference_genome.hpp"

#include <algorithm>
#include <iterator>
#include <utility>
#include <numeric>
#include <regex>
#include <stdexcept>

#include "fasta.hpp"
#include "caching_fasta.hpp"
#include "threadsafe_fasta.hpp"

ReferenceGenome::ReferenceGenome(std::unique_ptr<ReferenceGenomeImpl> impl)
:
impl_ {std::move(impl)},
name_ {impl_->get_reference_name()},
contig_names_(impl_->get_contig_names()),
contig_sizes_ {}
{
    for (const auto& contig_name : contig_names_) {
        contig_sizes_.emplace(contig_name, impl_->get_contig_size(contig_name));
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

bool ReferenceGenome::contains_region(const GenomicRegion& region) const noexcept
{
    return has_contig(region.get_contig_name()) && region.get_end() <= get_contig_size(region.get_contig_name());
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
    throw std::runtime_error {"get_contig_size: contig \"" + contig_name +
            "\" is not in reference genome \"" + name_ + "\""};
}

ReferenceGenome::SizeType ReferenceGenome::get_contig_size(const GenomicRegion& region) const
{
    return get_contig_size(region.get_contig_name());
}

GenomicRegion ReferenceGenome::get_contig_region(const std::string& contig_name) const
{
    return GenomicRegion {contig_name, 0, get_contig_size(contig_name)};
}

ReferenceGenome::SequenceType ReferenceGenome::get_sequence(const GenomicRegion& region) const
{
    return impl_->fetch_sequence(region);
}

// non-member functions

ReferenceGenome make_reference(boost::filesystem::path file_path, const std::size_t max_base_pair_cache,
                               const bool is_threaded)
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
    
    std::sort(std::begin(result), std::end(result),
              [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
    
    return result;
}

GenomicRegion::SizeType calculate_genome_size(const ReferenceGenome& reference)
{
    const auto contigs = reference.get_contig_names();
    return std::accumulate(std::cbegin(contigs), std::cend(contigs), GenomicRegion::SizeType {0},
                           [&reference] (const auto curr, const auto& contig) {
                               return curr + reference.get_contig_size(contig);
                           });
}

// Requires reference access to get contig sizes for partially specified regions (e.g. "4")
GenomicRegion parse_region(std::string region, const ReferenceGenome& reference)
{
    region.erase(std::remove(std::begin(region), std::end(region), ','), std::end(region));
    
    const static std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};
    
    std::smatch match;
    
    if (std::regex_search(region, match, re) && match.size() == 5) {
        auto contig_name = match.str(1);
        
        GenomicRegion::SizeType begin {}, end {};
        
        const auto contig_size = reference.get_contig_size(contig_name);
        
        if (match.str(2).empty()) {
            end = contig_size;
        } else {
            begin = static_cast<GenomicRegion::SizeType>(std::stoul(match.str(2)));
            
            if (match.str(3).empty()) {
                end = begin;
            } else if (match.str(4).empty()) {
                end = contig_size;
            } else {
                end = static_cast<GenomicRegion::SizeType>(std::stoul(match.str(4)));
            }
            
            if (begin > contig_size || end > contig_size) {
                throw std::runtime_error {"parse_region: given region (" + region +
                    ") larger than contig in " + reference.get_name()};
            }
        }
        
        return GenomicRegion {std::move(contig_name), begin, end};
    }
    
    throw std::runtime_error {"parse_region: given region (" + region + ") with invalid format"};
}
