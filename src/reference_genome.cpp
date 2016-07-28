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
#include "threadsafe_fasta.hpp"
#include "caching_fasta.hpp"

ReferenceGenome::ReferenceGenome(std::unique_ptr<ReferenceGenomeImpl> impl)
:
impl_ {std::move(impl)},
name_{},
contig_names_ {},
contig_sizes_ {}
{
    if (impl_->is_open()) {
        try {
            name_         = impl_->fetch_reference_name();
            contig_names_ = impl_->fetch_contig_names();
            
            contig_sizes_.reserve(contig_names_.size());
            
            for (const auto& contig_name : contig_names_) {
                contig_sizes_.emplace(contig_name, impl_->fetch_contig_size(contig_name));
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

const std::string& ReferenceGenome::name() const
{
    return name_;
}

bool ReferenceGenome::has_contig(const ContigName& contig) const noexcept
{
    return std::find(std::cbegin(contig_names_), std::cend(contig_names_), contig) != std::cend(contig_names_);
}

bool ReferenceGenome::contains_region(const GenomicRegion& region) const noexcept
{
    return this->has_contig(region.contig_name()) && region.end() <= this->contig_size(region.contig_name());
}

std::size_t ReferenceGenome::num_contigs() const noexcept
{
    return contig_names_.size();
}

const std::vector<ReferenceGenome::ContigName>& ReferenceGenome::contig_names() const noexcept
{
    return contig_names_;
}

ContigRegion::Size ReferenceGenome::contig_size(const ContigName& contig) const
{
    return contig_sizes_.at(contig);
}

ContigRegion::Size ReferenceGenome::contig_size(const GenomicRegion& region) const
{
    return this->contig_size(region.contig_name());
}

GenomicRegion ReferenceGenome::contig_region(const ContigName& contig) const
{
    return GenomicRegion {contig, GenomicRegion::Position {0}, this->contig_size(contig)};
}

ReferenceGenome::GeneticSequence ReferenceGenome::fetch_sequence(const GenomicRegion& region) const
{
    return impl_->fetch_sequence(region);
}

// non-member functions

ReferenceGenome make_reference(boost::filesystem::path reference_path,
                               const std::size_t max_cached_bases,
                               const bool is_threaded)
{
    std::unique_ptr<ReferenceGenomeImpl> impl_ {};
    
    try {
        if (is_threaded) {
            impl_ = std::make_unique<ThreadsafeFasta>(std::make_unique<Fasta>(reference_path));
        } else {
            impl_ = std::make_unique<Fasta>(std::move(reference_path));
        }
        
        if (max_cached_bases > 0) {
            double locality_bias {0.99}, forward_bias {0.99};
            
            if (is_threaded) {
                locality_bias = 0.25;
            }
            
            return ReferenceGenome {std::make_unique<CachingFasta>(std::move(impl_), max_cached_bases,
                                                                   locality_bias, forward_bias)};
        } else {
            return ReferenceGenome {std::move(impl_)};
        }
    } catch (const Fasta::MissingFastaIndex& e) {
        throw; // TODO: we could optionally make the index
    }
}

std::vector<GenomicRegion> get_all_contig_regions(const ReferenceGenome& reference)
{
    std::vector<GenomicRegion> result {};
    result.reserve(reference.num_contigs());
    
    for (const auto& contig : reference.contig_names()) {
        result.emplace_back(reference.contig_region(contig));
    }
    
    std::sort(std::begin(result), std::end(result),
              [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
    
    return result;
}

GenomicRegion::Position calculate_genome_size(const ReferenceGenome& reference)
{
    const auto& contigs = reference.contig_names();
    return std::accumulate(std::cbegin(contigs), std::cend(contigs),
                           GenomicRegion::Position {0},
                           [&reference] (const auto curr, const auto& contig) {
                               return curr + reference.contig_size(contig);
                           });
}

GenomicRegion parse_region(std::string region, const ReferenceGenome& reference)
{
    using Position = GenomicRegion::Position;
    
    region.erase(std::remove(std::begin(region), std::end(region), ','), std::end(region));
    
    static const std::regex re {"([^:]+)(?::(\\d+)(-)?(\\d*))?"};
    
    std::smatch match;
    
    if (std::regex_match(region, match, re) && match.size() == 5) {
        GenomicRegion::ContigName contig {match.str(1)};
        
        const auto contig_size = reference.contig_size(contig);
        
        Position begin {0}, end {0};
        
        if (match.str(2).empty()) {
            end = contig_size;
        } else {
            begin = static_cast<Position>(std::stoul(match.str(2)));
            
            if (match.str(3).empty()) {
                end = begin + 1;
            } else if (match.str(4).empty()) {
                end = contig_size;
            } else {
                end = static_cast<Position>(std::stoul(match.str(4)));
            }
            
            if (begin > end) {
                throw std::invalid_argument {"parse_region: given region ("
                        + region + ") with invalid format"};
            }
            
            if (begin > contig_size) {
                begin = contig_size;
            }
            
            if (end > contig_size) end = contig_size;
        }
        
        return GenomicRegion {std::move(contig), begin, end};
    }
    
    throw std::invalid_argument {"parse_region: given region (" + region + ") with invalid format"};
}
