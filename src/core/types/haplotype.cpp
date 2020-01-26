// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "haplotype.hpp"

#include <algorithm>
#include <iterator>
#include <stdexcept>
#include <iostream>
#include <cassert>

#include "io/reference/reference_genome.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable_range.hpp"
#include "utils/mappable_algorithms.hpp"
#include "variant.hpp"

namespace octopus {

template <typename T, typename M>
auto haplotype_overlap_range(const T& alleles, const M& mappable)
{
    return bases(overlap_range(alleles, mappable, BidirectionallySortedTag {}));
}

template <typename T, typename M>
auto haplotype_contained_range(const T& alleles, const M& mappable)
{
    return bases(contained_range(alleles, mappable));
}

// public methods

const GenomicRegion& Haplotype::mapped_region() const
{
    return region_;
}
    
namespace {

template <typename BidirIt, typename T>
BidirIt binary_find(BidirIt first, BidirIt last, const T& value)
{
    const auto itr = std::lower_bound(first, last, value);
    return (itr != last && *itr == value) ? itr : last;
}

} // namespace

bool Haplotype::contains(const ContigAllele& allele) const
{
    using octopus::contains; using std::cbegin; using std::cend;
    if (contains(region_.contig_region(), allele)) {
        if (begins_before(allele, explicit_allele_region_)) {
            if (is_before(allele, explicit_allele_region_)) {
                return allele.sequence() == fetch_reference_sequence(contig_region(allele));
            }
            const auto flank_region = left_overhang_region(explicit_allele_region_, contig_region(allele));
            if (copy(allele, flank_region).sequence() != fetch_reference_sequence(flank_region)) {
                return false;
            }
        }
        if (ends_before(explicit_allele_region_, allele)) {
            if (is_after(allele, explicit_allele_region_)) {
                return allele.sequence() == fetch_reference_sequence(contig_region(allele));
            }
            const auto flank_region = right_overhang_region(contig_region(allele), explicit_allele_region_);
            if (copy(allele, flank_region).sequence() != fetch_reference_sequence(flank_region)) {
                return false;
            }
        }
        const auto match_itr = binary_find(cbegin(explicit_alleles_), cend(explicit_alleles_), allele.mapped_region());
        if (match_itr != cend(explicit_alleles_)) {
            if (*match_itr == allele) return true;
            if (is_same_region(*match_itr, allele)) {
                // If the allele is not explcitly contained but the region is then it must be a different
                // allele, unless it is an insertion, in which case we must check the sequence
                if (is_insertion(allele)) {
                    return contains(*std::lower_bound(cbegin(explicit_alleles_), cend(explicit_alleles_),
                                                      allele.mapped_region()), allele);
                }
                return false;
            }
        }
        const auto overlapped = haplotype_overlap_range(explicit_alleles_, allele);
        if (overlapped.size() == 1 && contains(overlapped.front(), allele)) {
            return allele == copy(overlapped.front(), contig_region(allele));
        }
        return sequence(allele.mapped_region()) == allele.sequence();
    }
    return false;
}

bool Haplotype::contains(const Allele& allele) const
{
    if (!is_same_contig(allele, region_)) return false;
    return contains(demote(allele));
}

bool Haplotype::includes(const ContigAllele& allele) const
{
    using octopus::contains;
    if (!contains(region_.contig_region(), allele)) {
        return false;
    } else if (!explicit_alleles_.empty()) {
        if (contains(explicit_allele_region_, allele)) {
            return std::binary_search(std::cbegin(explicit_alleles_), std::cend(explicit_alleles_), allele);
        } else if (overlaps(explicit_allele_region_, allele)) {
            return false;
        } else if (is_after(allele, explicit_allele_region_)) {
            if (is_indel(allele)) return false;
            const auto ref_ritr = std::next(std::crbegin(sequence_), end_distance(allele, region_.contig_region()));
            assert(static_cast<std::size_t>(std::distance(ref_ritr, std::crend(sequence_))) >= allele.sequence().size());
            return std::equal(std::crbegin(allele.sequence()), std::crend(allele.sequence()), ref_ritr);
        }
    }
    if (is_indel(allele)) return false;
    const auto ref_itr = std::next(std::cbegin(sequence_), begin_distance(region_.contig_region(), allele));
    assert(static_cast<std::size_t>(std::distance(ref_itr, std::cend(sequence_))) >= allele.sequence().size());
    return std::equal(std::cbegin(allele.sequence()), std::cend(allele.sequence()), ref_itr);
}

bool Haplotype::includes(const Allele& allele) const
{
    if (!is_same_contig(allele, region_)) return false;
    return includes(demote(allele));
}

bool is_in_reference_flank(const ContigRegion& region, const ContigRegion& explicit_allele_region_,
                           const std::vector<ContigAllele>& explicit_alleles)
{
    if (overlaps(region, explicit_allele_region_)) {
        return false;
    }
    if (!are_adjacent(region, explicit_allele_region_)) {
        return true;
    }
    if (begins_before(region, explicit_allele_region_)) {
        return !is_simple_insertion(explicit_alleles.front());
    }
    return !is_simple_insertion(explicit_alleles.back());
}

Haplotype::NucleotideSequence Haplotype::sequence(const ContigRegion& region) const
{
    using std::cbegin; using std::cend; using octopus::contains;
    if (!contains(region_.contig_region(), region)) {
        throw std::out_of_range {"Haplotype: attempting to sequence from region not contained by Haplotype region"};
    }
    if (explicit_alleles_.empty()) {
        return sequence_.substr(begin_distance(region_.contig_region(), region), region_size(region));
    }
    if (is_in_reference_flank(region, explicit_allele_region_, explicit_alleles_)) {
        return fetch_reference_sequence(region);
    }
    NucleotideSequence result {};
    result.reserve(region_size(region)); // may be more or less depending on indels
    if (begins_before(region, explicit_allele_region_)) {
        append_reference(result, left_overhang_region(region, explicit_allele_region_));
    }
    auto overlapped_explicit_alleles = haplotype_overlap_range(explicit_alleles_, region);
    assert(!overlapped_explicit_alleles.empty());
    if (contains(overlapped_explicit_alleles.front(), region)) {
        append(result, copy(overlapped_explicit_alleles.front(), region));
        overlapped_explicit_alleles.advance_begin(1);
        if (!overlapped_explicit_alleles.empty() && is_empty_region(overlapped_explicit_alleles.front())) {
            append(result, overlapped_explicit_alleles.front());
        }
        return result;
    } else if (begins_before(overlapped_explicit_alleles.front(), region)) {
        append(result, copy(overlapped_explicit_alleles.front(),
                            *overlapped_region(overlapped_explicit_alleles.front(), region)));
        overlapped_explicit_alleles.advance_begin(1);
        if (overlapped_explicit_alleles.empty()) {
            append_reference(result, right_overhang_region(region, explicit_allele_region_));
            return result;
        }
    }
    bool region_ends_before_last_overlapped_allele {ends_before(region, overlapped_explicit_alleles.back())};
    if (region_ends_before_last_overlapped_allele) {
        overlapped_explicit_alleles.advance_end(-1); // as we don't want all of the last allele
        region_ends_before_last_overlapped_allele = true;
    }
    append(result, cbegin(overlapped_explicit_alleles), cend(overlapped_explicit_alleles));
    if (region_ends_before_last_overlapped_allele) {
        overlapped_explicit_alleles.advance_end(1); // as we previously removed this allele
        append(result, copy(overlapped_explicit_alleles.back(),
                            *overlapped_region(overlapped_explicit_alleles.back(), region)));
    } else if (ends_before(explicit_allele_region_, region)) {
        append_reference(result, right_overhang_region(region, explicit_allele_region_));
    }
    return result;
}

Haplotype::NucleotideSequence Haplotype::sequence(const GenomicRegion& region) const
{
    if (!is_same_contig(region, region_)) {
        throw std::logic_error {"Haplotype: cannot sequence from different contig"};
    }
    return sequence(region.contig_region());
}

const Haplotype::NucleotideSequence& Haplotype::sequence() const noexcept
{
    return sequence_;
}

Haplotype::NucleotideSequence::size_type Haplotype::sequence_size(const ContigRegion& region) const
{
    return sequence(region).size(); // TODO: can be improved
}

Haplotype::NucleotideSequence::size_type Haplotype::sequence_size(const GenomicRegion& region) const
{
    if (!is_same_contig(region, region_)) return 0;
    return sequence_size(region.contig_region());
}

std::vector<Variant> Haplotype::difference(const Haplotype& other) const
{
    std::vector<Variant> result {};
    result.reserve(explicit_alleles_.size());
    const auto& contig = region_.contig_name();
    for (const auto& allele : explicit_alleles_) {
        if (!other.contains(allele)) {
            result.emplace_back(GenomicRegion {contig, allele.mapped_region()},
                                other.sequence(allele.mapped_region()),
                                allele.sequence());
        }
    }
    return result;
}

CigarString Haplotype::cigar() const
{
    using Flag = CigarOperation::Flag;
    CigarString result {};
    if (!explicit_alleles_.empty()) {
        const auto reference = reference_.get().fetch_sequence(GenomicRegion {region_.contig_name(), explicit_allele_region_});
        result.reserve(2 * explicit_alleles_.size() + 2);
        auto curr_op_size = begin_distance(region_.contig_region(), explicit_allele_region_);
        auto curr_op_flag = Flag::sequenceMatch;
        for (std::size_t i {0}; i < explicit_alleles_.size(); ++i) {
            const auto& allele = explicit_alleles_[i];
            if (i > 0) {
                const auto& prev_allele = explicit_alleles_[i - 1];
                if (!are_adjacent(prev_allele, allele)) {
                    result.emplace_back(curr_op_size, curr_op_flag);
                    curr_op_flag = Flag::sequenceMatch;
                    curr_op_size = intervening_region_size(prev_allele, allele);
                }
            }
            auto allele_op_flag = curr_op_flag;
            CigarOperation::Size allele_op_size {0};
            if (is_indel(allele)) {
                if (is_simple_insertion(allele)) {
                    allele_op_flag = Flag::insertion;
                    allele_op_size += allele.sequence().size();
                } else if (is_simple_deletion(allele)) {
                    allele_op_flag = Flag::deletion;
                    allele_op_size += region_size(allele);
                } else {
                    // all complex indels are treated as replacements
                    if (curr_op_flag == Flag::deletion) {
                        curr_op_size += region_size(allele);
                        result.emplace_back(curr_op_size, curr_op_flag);
                        curr_op_flag = Flag::insertion;
                        curr_op_size = allele.sequence().size();
                    } else if (curr_op_flag == Flag::insertion) {
                        curr_op_size += allele.sequence().size();
                        result.emplace_back(curr_op_size, curr_op_flag);
                        curr_op_flag = Flag::deletion;
                        curr_op_size = region_size(allele);
                    } else {
                        if (curr_op_size > 0) {
                            result.emplace_back(curr_op_size, curr_op_flag);
                        }
                        result.emplace_back(region_size(allele), Flag::deletion);
                        curr_op_flag = Flag::insertion;
                        curr_op_size = allele.sequence().size();
                    }
                }
            } else if (!is_empty_region(allele)) {
                const auto ref_idx = static_cast<std::size_t>(begin_distance(explicit_allele_region_, allele));
                assert(ref_idx < reference.size());
                if (region_size(allele) == 1) {
                    if (allele.sequence()[0] == reference[ref_idx]) {
                        allele_op_flag = Flag::sequenceMatch;
                    } else {
                        allele_op_flag = Flag::substitution;
                    }
                    ++allele_op_size;
                } else {
                    assert(reference.size() >= ref_idx + allele.sequence().size());
                    if (std::equal(std::cbegin(allele.sequence()), std::cend(allele.sequence()),
                                   std::next(std::cbegin(reference), ref_idx))) {
                        allele_op_flag = Flag::sequenceMatch;
                    } else {
                        allele_op_flag = Flag::alignmentMatch;
                    }
                    allele_op_size += region_size(allele);
                }
            }
            if (allele_op_flag == curr_op_flag) {
                curr_op_size += allele_op_size;
            } else {
                if (curr_op_size > 0) {
                    result.emplace_back(curr_op_size, curr_op_flag);
                }
                curr_op_flag = allele_op_flag;
                curr_op_size = allele_op_size;
            }
        }
        const auto rhs_ref_flank_size = end_distance(explicit_allele_region_, region_.contig_region());
        if (curr_op_flag == Flag::sequenceMatch) {
            curr_op_size += rhs_ref_flank_size;
        } else {
            result.emplace_back(curr_op_size, curr_op_flag);
            curr_op_flag = Flag::sequenceMatch;
            curr_op_size = rhs_ref_flank_size;
        }
        if (curr_op_size > 0) {
            result.emplace_back(curr_op_size, curr_op_flag);
        }
    } else {
        result.emplace_back(size(region_), Flag::sequenceMatch);
    }
    assert(octopus::sequence_size(result) == sequence_.size());
    assert(reference_size(result) == size(region_));
    return result;
}

std::size_t Haplotype::get_hash() const noexcept
{
    return cached_hash_;
}

std::pair<Haplotype::AlleleIterator, Haplotype::AlleleIterator>
Haplotype::alleles() const noexcept
{
    return std::make_pair(std::cbegin(explicit_alleles_), std::cend(explicit_alleles_));
}

// private methods

void Haplotype::append(NucleotideSequence& result, const ContigAllele& allele) const
{
    result.append(allele.sequence());
}

void Haplotype::append(NucleotideSequence& result, AlleleIterator first, AlleleIterator last) const
{
    std::for_each(first, last, [this, &result] (const ContigAllele& allele) { append(result, allele); });
}

void Haplotype::append_reference(NucleotideSequence& result, const ContigRegion& region) const
{
    if (is_before(region, explicit_allele_region_)) {
        const auto offset = begin_distance(region_.contig_region(), region);
        const auto it = std::next(std::cbegin(sequence_), offset);
        result.append(it, std::next(it, region_size(region)));
    } else {
        const auto offset = end_distance(region, region_.contig_region());
        const auto it = std::prev(std::cend(sequence_), offset);
        result.append(std::prev(it, region_size(region)), it);
    }
}

Haplotype::NucleotideSequence Haplotype::fetch_reference_sequence(const ContigRegion& region) const
{
    NucleotideSequence result {};
    result.reserve(region_size(region));
    append_reference(result, region);
    return result;
}

// Builder

Haplotype::Builder::Builder(const GenomicRegion& region, const ReferenceGenome& reference)
:
region_ {region},
reference_ {reference}
{}

bool Haplotype::Builder::can_push_back(const ContigAllele& allele) const noexcept
{
    return explicit_alleles_.empty() || is_after(allele, explicit_alleles_.back());
}

bool Haplotype::Builder::can_push_back(const Allele& allele) const noexcept
{
    return is_same_contig(allele, region_) && (explicit_alleles_.empty() || is_after(contig_region(allele), explicit_alleles_.back()));
}

bool Haplotype::Builder::can_push_front(const ContigAllele& allele) const noexcept
{
    return explicit_alleles_.empty() || is_after(explicit_alleles_.front(), allele);
}

bool Haplotype::Builder::can_push_front(const Allele& allele) const noexcept
{
    return is_same_contig(allele, region_) && (explicit_alleles_.empty() || is_after(explicit_alleles_.front(), contig_region(allele)));
}

void Haplotype::Builder::push_back(const ContigAllele& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(allele, explicit_alleles_.back())) {
            throw std::logic_error {"Haplotype::push_back called with out-of-order Allele"};
        }
        if (!are_adjacent(explicit_alleles_.back(), allele)) {
            explicit_alleles_.push_back(get_intervening_reference_allele(explicit_alleles_.back(), allele));
        }
    }
    update_region(allele);
    explicit_alleles_.push_back(allele);
}

void Haplotype::Builder::push_back(ContigAllele&& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(allele, explicit_alleles_.back())) {
            throw std::logic_error {"Haplotype::push_back called with out-of-order Allele"};
        }
        if (!are_adjacent(explicit_alleles_.back(), allele)) {
            explicit_alleles_.push_back(get_intervening_reference_allele(explicit_alleles_.back(), allele));
        }
    }
    update_region(allele);
    explicit_alleles_.push_back(std::move(allele));
}

void Haplotype::Builder::push_front(const ContigAllele& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(explicit_alleles_.front(), allele)) {
            throw std::logic_error {"Haplotype::push_front called with out-of-order Allele"};
        }
        if (!are_adjacent(allele, explicit_alleles_.front())) {
            explicit_alleles_.push_front(get_intervening_reference_allele(allele, explicit_alleles_.front()));
        }
    }
    update_region(allele);
    explicit_alleles_.push_front(allele);
}

void Haplotype::Builder::push_front(ContigAllele&& allele)
{
    if (!explicit_alleles_.empty()) {
        if (!is_after(explicit_alleles_.front(), allele)) {
            throw std::logic_error {"Haplotype::push_front called with out-of-order Allele"};
        }
        if (!are_adjacent(allele, explicit_alleles_.front())) {
            explicit_alleles_.push_front(get_intervening_reference_allele(allele, explicit_alleles_.front()));
        }
    }
    update_region(allele);
    explicit_alleles_.push_front(std::move(allele));
}

void Haplotype::Builder::push_back(const Allele& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::domain_error {"Haplotype::push_back called with Allele on different contig"};
    }
    push_back(ContigAllele {contig_region(allele), allele.sequence()});
}

void Haplotype::Builder::push_front(const Allele& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::domain_error {"Haplotype::push_front called with Allele on different contig"};
    }
    push_front(ContigAllele {contig_region(allele), allele.sequence()});
}

void Haplotype::Builder::push_back(Allele&& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::domain_error {"Haplotype::push_back called with Allele on different contig"};
    }
    push_back(ContigAllele {contig_region(allele), allele.sequence()});
}

void Haplotype::Builder::push_front(Allele&& allele)
{
    if (!is_same_contig(allele, region_)) {
        throw std::domain_error {"Haplotype::push_front called with Allele on different contig"};
    }
    push_front(ContigAllele {contig_region(allele), allele.sequence()});
}

Haplotype Haplotype::Builder::build()
{
    return Haplotype {
        std::move(region_),
        std::make_move_iterator(std::begin(explicit_alleles_)),
        std::make_move_iterator(std::end(explicit_alleles_)),
        reference_
    };
}

void Haplotype::Builder::update_region(const ContigAllele& allele) noexcept
{
    const auto new_contig_region = encompassing_region(region_.contig_region(), allele);
    region_ = GenomicRegion {region_.contig_name(), new_contig_region};
}

void Haplotype::Builder::update_region(const Allele& allele)
{
    region_ = encompassing_region(region_, allele);
}

ContigAllele Haplotype::Builder::get_intervening_reference_allele(const ContigAllele& lhs, const ContigAllele& rhs) const
{
    const auto region = *intervening_region(lhs, rhs);
    return ContigAllele {region, reference_.get().fetch_sequence(GenomicRegion {region_.contig_name(), region})};
}

// non-member methods

Haplotype::NucleotideSequence::size_type sequence_size(const Haplotype& haplotype) noexcept
{
    return haplotype.sequence().size();
}

bool is_sequence_empty(const Haplotype& haplotype) noexcept
{
    return haplotype.sequence().empty();
}

bool contains(const Haplotype& lhs, const Allele& rhs)
{
    return lhs.contains(rhs);
}

bool contains(const Haplotype& lhs, const Haplotype& rhs)
{
    return contains(mapped_region(lhs), mapped_region(rhs)) && lhs.sequence(rhs.region_) == rhs.sequence();
}

bool includes(const Haplotype& lhs, const Allele& rhs)
{
    return lhs.includes(rhs);
}

namespace detail {

Haplotype do_copy(const Haplotype& haplotype, const GenomicRegion& region, std::true_type)
{
    using std::end; using std::cbegin; using std::cend; using std::prev;
    if (!is_same_contig(haplotype, region)) {
        throw std::logic_error {"Haplotype: trying to copy region from different contig"};
    }
    if (!contains(contig_region(haplotype), contig_region(region))) {
        throw std::logic_error {"Haplotype: trying to copy uncontained region"};
    }
    if (is_same_region(haplotype, region)) return haplotype;
    Haplotype::Builder result {region, haplotype.reference_};
    if (haplotype.explicit_alleles_.empty()) return result.build();
    const auto& contig_region = region.contig_region();
    if (contains(contig_region, haplotype.explicit_allele_region_)) {
        result.explicit_alleles_.assign(cbegin(haplotype.explicit_alleles_), cend(haplotype.explicit_alleles_));
        return result.build();
    }
    if (!overlaps(contig_region, haplotype.explicit_allele_region_)) return result.build();
    auto overlapped = haplotype_overlap_range(haplotype.explicit_alleles_, region.contig_region());
    assert(!overlapped.empty());
    if (is_empty(contig_region)) {
        if (!is_empty_region(overlapped.front()) && are_adjacent(contig_region, overlapped.front())) {
            overlapped.advance_begin(1);
        }
        if (!overlapped.empty() && is_empty_region(overlapped.front())) {
            result.push_back(overlapped.front());
        } else {
            result.push_back(ContigAllele {contig_region, ""});
        }
        return result.build();
    }
    if (!contains(contig_region, overlapped.front())) {
        result.push_front(copy(overlapped.front(), *overlapped_region(overlapped.front(), contig_region)));
        overlapped.advance_begin(1);
    }
    if (!overlapped.empty()) {
        if (contains(contig_region, overlapped.back())) {
            result.explicit_alleles_.insert(end(result.explicit_alleles_), cbegin(overlapped), cend(overlapped));
        } else {
            result.explicit_alleles_.insert(end(result.explicit_alleles_), cbegin(overlapped), prev(cend(overlapped)));
            result.push_back(copy(overlapped.back(), *overlapped_region(overlapped.back(), contig_region)));
        }
    }
    return result.build();
}

Allele do_copy(const Haplotype& haplotype, const GenomicRegion& region, std::false_type)
{
    return Allele {region, haplotype.sequence(region)};
}

} // namespace detail

ContigAllele copy(const Haplotype& haplotype, const ContigRegion& region)
{
    return ContigAllele {region, haplotype.sequence(region)};
}

bool is_reference(const Haplotype& haplotype)
{
    if (haplotype.explicit_alleles_.empty()) return true;
    return haplotype.sequence() == haplotype.reference_.get().fetch_sequence(haplotype.mapped_region());
}

Haplotype expand(const Haplotype& haplotype, Haplotype::MappingDomain::Size n)
{
    if (n == 0) return haplotype;
    return Haplotype {
        expand(mapped_region(haplotype), n),
        std::cbegin(haplotype.explicit_alleles_), std::cend(haplotype.explicit_alleles_),
        haplotype.reference_
    };
}

Haplotype remap(const Haplotype& haplotype, const GenomicRegion& region)
{
    if (is_same_region(haplotype, region)) {
        return haplotype;
    } else if (contains(region, haplotype)) {
        return Haplotype {
            region, std::cbegin(haplotype.explicit_alleles_), std::cend(haplotype.explicit_alleles_), haplotype.reference_
        };
    } else if (contains(haplotype, region)) {
        return copy<Haplotype>(haplotype, region);
    } else if (is_same_contig(haplotype, region)) {
        const auto remap_alleles = haplotype_contained_range(haplotype.explicit_alleles_, region.contig_region());
        return Haplotype {
            region, std::cbegin(remap_alleles), std::cend(remap_alleles), haplotype.reference_
        };
    } else {
        return Haplotype {region, haplotype.reference_};
    }
}

std::vector<Variant> difference(const Haplotype& lhs, const Haplotype& rhs)
{
    auto result = lhs.difference(rhs);
    auto diffs2 = rhs.difference(lhs);
    const auto itr = result.insert(std::end(result),
                                   std::make_move_iterator(std::begin(diffs2)),
                                   std::make_move_iterator(std::end(diffs2)));
    std::inplace_merge(std::begin(result), itr, std::end(result));
    return result;
}

bool operator==(const Haplotype& lhs, const Haplotype& rhs)
{
    return lhs.mapped_region() == rhs.mapped_region() && lhs.sequence() == rhs.sequence();
}

bool operator<(const Haplotype& lhs, const Haplotype& rhs)
{
    return lhs.mapped_region() == rhs.mapped_region() ? lhs.sequence() < rhs.sequence() : lhs.mapped_region() < rhs.mapped_region();
}

bool HaveSameAlleles::operator()(const Haplotype &lhs, const Haplotype &rhs) const
{
    return lhs.explicit_alleles_ == rhs.explicit_alleles_;
}

bool have_same_alleles(const Haplotype& lhs, const Haplotype& rhs)
{
    return HaveSameAlleles()(lhs, rhs);
}

IsLessComplex::IsLessComplex(boost::optional<Haplotype> reference) : reference_ {std::move(reference)} {}

bool IsLessComplex::operator()(const Haplotype& lhs, const Haplotype& rhs) const
{
    if (lhs.explicit_alleles_.size() != rhs.explicit_alleles_.size()) {
        return lhs.explicit_alleles_.size() < rhs.explicit_alleles_.size();
    }
    if (reference_) {
        return lhs.difference(*reference_).size() < rhs.difference(*reference_).size();
    }
    // otherwise prefer the sequence with the least amount of indels
    auto score = std::inner_product(std::cbegin(lhs.explicit_alleles_), std::cend(lhs.explicit_alleles_),
                                    std::cbegin(rhs.explicit_alleles_), 0, std::plus<> {},
                                    [] (const auto& lhs, const auto& rhs) {
                                        if (lhs == rhs) {
                                            return 0;
                                        } else if (is_indel(lhs)) {
                                            if (is_indel(rhs)) {
                                                return 0;
                                            } else {
                                                return -1;
                                            }
                                        } else if (is_indel(rhs)) {
                                            return 1;
                                        } else {
                                            return 0;
                                        }
                                    });
    return score >= 0;
}

bool StrictLess::operator()(const Haplotype& lhs, const Haplotype& rhs) const
{
    if (lhs.mapped_region() == rhs.mapped_region()) {
        if (lhs.sequence() != rhs.sequence()) {
            return lhs.sequence() < rhs.sequence();
        } else {
            return lhs.explicit_alleles_ < rhs.explicit_alleles_;
        }
    } else {
        return lhs.mapped_region() < rhs.mapped_region();
    }
}

bool are_equal_in_region(const Haplotype& lhs, const Haplotype& rhs, const GenomicRegion& region)
{
    return copy<Haplotype>(lhs, region) == copy<Haplotype>(rhs, region);
}

std::ostream& operator<<(std::ostream& os, const Haplotype& haplotype)
{
    os << haplotype.mapped_region() << " " << haplotype.sequence();
    return os;
}

namespace debug {

void print_alleles(const Haplotype& haplotype)
{
    print_alleles(std::cout, haplotype);
}

void print_variant_alleles(const Haplotype& haplotype)
{
    print_variant_alleles(std::cout, haplotype);
}

} // namespace debug
} // namespace octopus
