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
#include "threadsafe_caching_fasta.hpp"

ReferenceGenome::ReferenceGenome(std::unique_ptr<ReferenceGenomeImpl> impl)
:
impl_ {std::move(impl)},
name_{},
contig_names_ {},
contig_sizes_ {}
{
    if (impl_->is_open()) {
        try {
            name_         = impl_->get_reference_name();
            contig_names_ = impl_->get_contig_names();
            
            contig_sizes_.reserve(contig_names_.size());
            
            for (const auto& contig_name : contig_names_) {
                contig_sizes_.emplace(contig_name, impl_->get_contig_size(contig_name));
            }
        } catch (...) {
            impl_.reset(nullptr);
        }
    } else {
        impl_.reset(nullptr); // no point in keeping bad handle
    }
}

bool ReferenceGenome::is_good() const noexcept
{
    return impl_ != nullptr;
}

const std::string& ReferenceGenome::get_name() const
{
    return name_;
}

bool ReferenceGenome::has_contig(const ContigNameType& contig) const noexcept
{
    return std::find(std::cbegin(contig_names_), std::cend(contig_names_), contig) != std::cend(contig_names_);
}

bool ReferenceGenome::contains_region(const GenomicRegion& region) const noexcept
{
    return has_contig(region.get_contig_name()) && region.get_end() <= get_contig_size(region.get_contig_name());
}

std::size_t ReferenceGenome::num_contigs() const noexcept
{
    return contig_names_.size();
}

const std::vector<ReferenceGenome::ContigNameType>& ReferenceGenome::get_contig_names() const noexcept
{
    return contig_names_;
}

ReferenceGenome::SizeType ReferenceGenome::get_contig_size(const ContigNameType& contig) const
{
    return contig_sizes_.at(contig);
}

ReferenceGenome::SizeType ReferenceGenome::get_contig_size(const GenomicRegion& region) const
{
    return get_contig_size(region.get_contig_name());
}

GenomicRegion ReferenceGenome::get_contig_region(const ContigNameType& contig) const
{
    return GenomicRegion {contig, GenomicRegion::SizeType {0}, get_contig_size(contig)};
}

ReferenceGenome::SequenceType ReferenceGenome::get_sequence(const GenomicRegion& region) const
{
    return impl_->fetch_sequence(region);
}

// non-member functions

ReferenceGenome make_reference(boost::filesystem::path reference_path,
                               const std::size_t max_base_pair_cache,
                               const bool is_threaded)
{
    if (is_threaded) {
        if (max_base_pair_cache > 0) {
            return ReferenceGenome {std::make_unique<ThreadsafeCachingFasta>(std::move(reference_path),
                                                                             max_base_pair_cache)};
        } else {
            return ReferenceGenome {std::make_unique<ThreadsafeFasta>(std::move(reference_path))};
        }
    } else {
        if (max_base_pair_cache > 0) {
            return ReferenceGenome {std::make_unique<CachingFasta>(std::move(reference_path),
                                                                   max_base_pair_cache)};
        } else {
            return ReferenceGenome {std::make_unique<Fasta>(std::move(reference_path))};
        }
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

boost::optional<GenomicRegion> parse_region(std::string region, const ReferenceGenome& reference)
{
    region.erase(std::remove(std::begin(region), std::end(region), ','), std::end(region));
    
    const static std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};
    
    std::smatch match;
    
    if (std::regex_search(region, match, re) && match.size() == 5) {
        GenomicRegion::ContigNameType contig {match.str(1)};
        
        if (!reference.has_contig(contig)) return boost::none;
        
        const auto contig_size = reference.get_contig_size(contig);
        
        GenomicRegion::SizeType begin {0}, end {0};
        
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
            
            if (begin > contig_size) return boost::none;
            
            if (end > contig_size) end = contig_size;
        }
        
        return GenomicRegion {std::move(contig), begin, end};
    }
    
    return boost::none;
}
