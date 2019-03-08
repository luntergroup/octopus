// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_filter_hpp
#define read_filter_hpp

#include <string>
#include <utility>
#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

#include "basics/cigar_string.hpp"
#include "basics/aligned_read.hpp"
#include "basics/mappable_reference_wrapper.hpp"

namespace octopus { namespace readpipe
{
// All filters are nameable

class Nameable
{
    std::string name_;
    
public:
    Nameable() = delete;
    
    Nameable(std::string name) : name_ {std::move(name)} {}
    
    const std::string& name() const noexcept
    {
        return name_;
    }
};

// Basic filters

class BasicReadFilter : public Nameable
{
public:
    BasicReadFilter() = delete;
    
    virtual ~BasicReadFilter() = default;
    
    bool operator()(const AlignedRead& read) const noexcept
    {
        return passes(read);
    }
    
protected:
    BasicReadFilter(std::string name) : Nameable {std::move(name)} {};
    
private:
    virtual bool passes(const AlignedRead&) const noexcept = 0;
};

struct HasWellFormedCigar : BasicReadFilter
{
    HasWellFormedCigar();
    HasWellFormedCigar(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct HasValidBaseQualities : BasicReadFilter
{
    HasValidBaseQualities();
    HasValidBaseQualities(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotSecondaryAlignment : BasicReadFilter
{
    IsNotSecondaryAlignment();
    IsNotSecondaryAlignment(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotSupplementaryAlignment : BasicReadFilter
{
    IsNotSupplementaryAlignment();
    IsNotSupplementaryAlignment(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsGoodMappingQuality : BasicReadFilter
{
    using MappingQuality = AlignedRead::MappingQuality;
    
    IsGoodMappingQuality() = delete;
    
    explicit IsGoodMappingQuality(MappingQuality good_mapping_quality);
    
    IsGoodMappingQuality(std::string name, MappingQuality good_mapping_quality);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    MappingQuality good_mapping_quality_;
};

struct HasSufficientGoodBaseFraction : BasicReadFilter
{
    using BaseQuality = AlignedRead::BaseQuality;
    
    HasSufficientGoodBaseFraction() = delete;
    explicit HasSufficientGoodBaseFraction(BaseQuality good_base_quality,
                                           double min_good_base_fraction);
    HasSufficientGoodBaseFraction(std::string name, BaseQuality good_base_quality,
                                  double min_good_base_fraction);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    BaseQuality good_base_quality_;
    double min_good_base_fraction_;
};

struct HasSufficientGoodQualityBases : BasicReadFilter
{
    using BaseQuality = AlignedRead::BaseQuality;
    
    HasSufficientGoodQualityBases() = delete;
    explicit HasSufficientGoodQualityBases(BaseQuality good_base_quality,
                                           unsigned min_good_bases);
    HasSufficientGoodQualityBases(std::string name, BaseQuality good_base_quality,
                                  unsigned min_good_bases);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    BaseQuality good_base_quality_;
    unsigned min_good_bases_;
};

struct IsMapped : BasicReadFilter
{
    IsMapped();
    IsMapped(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotChimeric : BasicReadFilter
{
    IsNotChimeric();
    IsNotChimeric(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNextSegmentMapped : BasicReadFilter
{
    IsNextSegmentMapped();
    IsNextSegmentMapped(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotMarkedDuplicate : BasicReadFilter
{
    IsNotMarkedDuplicate();
    IsNotMarkedDuplicate(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsShort : BasicReadFilter
{
    using Length = AlignedRead::NucleotideSequence::size_type;
    
    IsShort() = delete;
    explicit IsShort(Length max_length);
    IsShort(std::string name, Length max_length);
    
    bool passes(const AlignedRead& read) const noexcept override;

private:
    Length max_length_;
};

struct IsLong : BasicReadFilter
{
    using Length = AlignedRead::AlignedRead::NucleotideSequence::size_type;
    
    IsLong() = delete;
    explicit IsLong(Length min_length);
    IsLong(std::string name, Length min_length);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    Length min_length_;
};

struct IsNotContaminated : BasicReadFilter
{
    IsNotContaminated();
    IsNotContaminated(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotMarkedQcFail : BasicReadFilter
{
    IsNotMarkedQcFail();
    IsNotMarkedQcFail(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsProperTemplate : BasicReadFilter
{
    IsProperTemplate();
    IsProperTemplate(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};
    
struct IsLocalTemplate : BasicReadFilter
{
    IsLocalTemplate();
    IsLocalTemplate(std::string name);
    
    bool passes(const AlignedRead& read) const noexcept override;
};

// Context filters

template <typename BidirIt>
class ContextReadFilter : public Nameable
{
public:
    ContextReadFilter() = delete;
    
    virtual ~ContextReadFilter() = default;
    
    BidirIt remove(BidirIt first, BidirIt last) const
    {
        return do_remove(first, last);
    }
    
    BidirIt partition(BidirIt first, BidirIt last) const
    {
        return do_partition(first, last);
    }
    
protected:
    ContextReadFilter(std::string name) : Nameable {std::move(name)} {};
    
private:
    virtual BidirIt do_remove(BidirIt first, BidirIt last) const = 0;
    virtual BidirIt do_partition(BidirIt first, BidirIt last) const = 0;
};

bool primary_segments_are_duplicates(const AlignedRead& lhs, const AlignedRead& rhs) noexcept;
bool other_segments_are_duplicates(const AlignedRead& lhs, const AlignedRead& rhs) noexcept;

struct IsDuplicate
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept;
};

template <typename Range>
bool no_duplicates(const AlignedRead& read, const Range& reads) noexcept
{
    return std::none_of(std::cbegin(reads), std::cend(reads), [&read] (auto other_itr) { return IsDuplicate{}(read, *other_itr); });
}

template <typename ForwardIt>
struct IsNotDuplicate : ContextReadFilter<ForwardIt>
{
    IsNotDuplicate() : ContextReadFilter<ForwardIt> {"IsNotOctopusDuplicate"} {}
    IsNotDuplicate(std::string name)
    : ContextReadFilter<ForwardIt> {std::move(name)} {}
    
    ForwardIt do_remove(ForwardIt first, ForwardIt last) const override
    {
        // Recall that reads come sorted w.r.t operator< and it is therefore not guaranteed that 'duplicate'
        // reads (according to IsDuplicate) will be adjacent to one another. In particular, operator< only
        // guarantees that the read segment described in the A
        first = std::adjacent_find(first, last, [] (const auto& lhs, const auto& rhs) { return primary_segments_are_duplicates(lhs, rhs); });
        if (first != last) {
            std::vector<ForwardIt> buffer {first++};
            buffer.reserve(100);
            for (auto itr = first; itr != last; ++itr) {
                if (primary_segments_are_duplicates(*itr, *buffer.front())) { // can check any read in buffer
                    if (no_duplicates(*itr, buffer)) {
                        if (itr != first) *first = std::move(*itr);
                        buffer.emplace_back(first++);
                    }
                } else {
                    if (itr != first) *first = std::move(*itr);
                    buffer.assign({first++});
                }
            }
        }
        return first;
    }
    
    ForwardIt do_partition(ForwardIt first, ForwardIt last) const override
    {
        // TODO: we need a clever stable_partition_unique implementation.
        // See my question:
        // http://stackoverflow.com/questions/36888033/implementing-partition-unique-and-stable-partition-unique-algorithms
        return last;
    }
};

} // namespace readpipe
} // namespace octopus

#endif
