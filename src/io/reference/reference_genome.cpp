// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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

namespace octopus {

ReferenceGenome::ReferenceGenome(std::unique_ptr<io::ReferenceReader> impl)
:
impl_ {std::move(impl)},
name_{},
contig_sizes_ {}
{
    if (impl_->is_open()) {
        try {
            name_ = impl_->fetch_reference_name();
            
            auto contig_names = impl_->fetch_contig_names();
            
            contig_sizes_.reserve(contig_names.size());
            
            for (const auto& contig_name : contig_names) {
                contig_sizes_.emplace(contig_name, impl_->fetch_contig_size(contig_name));
            }
        } catch (...) {
            impl_.reset(nullptr);
        }
    } else {
        impl_.reset(nullptr); // no point in keeping bad handle
    }
}

ReferenceGenome::ReferenceGenome(const ReferenceGenome& other)
:
impl_ {other.impl_->clone()},
name_ {other.name_},
contig_sizes_ {other.contig_sizes_}
{}

ReferenceGenome& ReferenceGenome::operator=(ReferenceGenome other)
{
    using std::swap;
    swap(impl_, other.impl_);
    swap(name_, other.name_);
    swap(contig_sizes_, other.contig_sizes_);
    return *this;
}

const std::string& ReferenceGenome::name() const
{
    return name_;
}

bool ReferenceGenome::has_contig(const ContigName& contig) const noexcept
{
    return contig_sizes_.count(contig) == 1;
}

std::size_t ReferenceGenome::num_contigs() const noexcept
{
    return contig_sizes_.size();
}

std::vector<ReferenceGenome::ContigName> ReferenceGenome::contig_names() const
{
    std::vector<ContigName> result {};
    result.reserve(contig_sizes_.size());
    
    std::transform(std::cbegin(contig_sizes_), std::cend(contig_sizes_), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    
    return result;
}

ContigRegion::Size ReferenceGenome::contig_size(const ContigName& contig) const
{
    return contig_sizes_.at(contig);
}

GenomicRegion ReferenceGenome::contig_region(const ContigName& contig) const
{
    return GenomicRegion {contig, GenomicRegion::Position {0}, this->contig_size(contig)};
}

bool ReferenceGenome::contains(const GenomicRegion& region) const noexcept
{
    return this->has_contig(region.contig_name()) && region.end() <= this->contig_size(region.contig_name());
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
    using namespace io;
    
    std::unique_ptr<ReferenceReader> impl_ {};
    
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
        
        if (match.length(2) == 0) {
            end = contig_size;
        } else {
            begin = static_cast<Position>(std::stoul(match.str(2)));
            
            if (match.length(3) == 0) {
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
    
} // namespace octopus
