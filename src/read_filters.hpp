//
//  read_filters.h
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filters__
#define __Octopus__read_filters__

#include <string>
#include <utility>
#include <algorithm>
#include <iterator>

#include "aligned_read.hpp"
#include "cigar_string.hpp"

namespace Octopus { namespace ReadFilters {

class BasicReadFilter
{
public:
    BasicReadFilter() = delete;
    
    virtual ~BasicReadFilter() = default;
    
    bool operator()(const AlignedRead& read) const noexcept
    {
        return passes(read);
    }
    
    const std::string& name() const noexcept
    {
        return name_;
    }
    
protected:
    BasicReadFilter(std::string name) : name_ {std::move(name)} {};
    
    std::string name_;
    
private:
    virtual bool passes(const AlignedRead&) const noexcept = 0;
};

struct IsNotSecondaryAlignment : BasicReadFilter
{
    IsNotSecondaryAlignment() : BasicReadFilter {"IsNotSecondaryAlignment"} {}
    
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotSupplementaryAlignment : BasicReadFilter
{
    IsNotSupplementaryAlignment() : BasicReadFilter {"IsNotSupplementaryAlignment"} {}
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsGoodMappingQuality : BasicReadFilter
{
    using QualityType = AlignedRead::QualityType;
    
    IsGoodMappingQuality() = delete;
    
    explicit IsGoodMappingQuality(QualityType good_mapping_quality);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    QualityType good_mapping_quality_;
};

struct HasSufficientGoodBaseFraction : BasicReadFilter
{
    using QualityType = AlignedRead::QualityType;
    
    HasSufficientGoodBaseFraction() = delete;
    
    explicit HasSufficientGoodBaseFraction(QualityType good_base_quality,
                                           double min_good_base_fraction);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    QualityType good_base_quality_;
    double min_good_base_fraction_;
};

struct HasSufficientGoodQualityBases : BasicReadFilter
{
    using QualityType = AlignedRead::QualityType;
    
    HasSufficientGoodQualityBases() = delete;
    
    explicit HasSufficientGoodQualityBases(QualityType good_base_quality,
                                           unsigned min_good_bases);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    QualityType good_base_quality_;
    unsigned min_good_bases_;
};

struct IsMapped : BasicReadFilter
{
    IsMapped() : BasicReadFilter {"IsMapped"} {};
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotChimeric : BasicReadFilter
{
    IsNotChimeric() : BasicReadFilter {"IsNotChimeric"} {}
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNextSegmentMapped : BasicReadFilter
{
    IsNextSegmentMapped() : BasicReadFilter {"IsNextSegmentMapped"} {}
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotMarkedDuplicate : BasicReadFilter
{
    IsNotMarkedDuplicate() : BasicReadFilter {"IsNotMarkedDuplicate"} {}
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsShort : BasicReadFilter
{
    using SizeType = AlignedRead::SizeType;
    
    IsShort() = delete;
    
    explicit IsShort(SizeType max_length);
    
    bool passes(const AlignedRead& read) const noexcept override;

private:
    SizeType max_length_;
};

struct IsLong : BasicReadFilter
{
    using SizeType = AlignedRead::SizeType;
    
    IsLong() = delete;
    
    explicit IsLong(SizeType max_length);
    
    bool passes(const AlignedRead& read) const noexcept override;
    
private:
    SizeType max_length_;
};

struct IsNotContaminated : BasicReadFilter
{
    IsNotContaminated() : BasicReadFilter {"IsNotContaminated"} {}
    bool passes(const AlignedRead& read) const noexcept override;
};

struct IsNotMarkedQcFail : BasicReadFilter
{
    IsNotMarkedQcFail() : BasicReadFilter {"IsNotMarkedQcFail"} {}
    bool passes(const AlignedRead& read) const noexcept override;
};

// Context-based filters

struct FilterDuplicates
{
    template <typename ForwardIt>
    ForwardIt operator()(ForwardIt first_read, ForwardIt last_read) const
    {
        return std::unique(first_read, last_read, IsDuplicate {});
    }
};

} // namespace ReadFilters
} // namespace Octopus

#endif /* defined(__Octopus__read_filters__) */
