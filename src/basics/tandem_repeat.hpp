// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef tandem_repeat_hpp
#define tandem_repeat_hpp

#include <string>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "basics/genomic_region.hpp"

namespace octopus {

class TandemRepeat : public Mappable<TandemRepeat>, public Comparable<TandemRepeat>
{
public:
    using NucleotideSequence = std::string;
    using SizeType = GenomicRegion::Size;
    
    TandemRepeat() = delete;
    template <typename R, typename S>
    TandemRepeat(R&& region, S&& motif)
    : region_ {std::forward<R>(region)}
    , motif_ {std::forward<S>(motif)}
    {}
    
    TandemRepeat(const TandemRepeat&)            = default;
    TandemRepeat& operator=(const TandemRepeat&) = default;
    TandemRepeat(TandemRepeat&&)                 = default;
    TandemRepeat& operator=(TandemRepeat&&)      = default;
    
    ~TandemRepeat() = default;
    
    SizeType period() const noexcept { return motif_.size(); }
    const NucleotideSequence& motif() const noexcept { return motif_; }
    const GenomicRegion& mapped_region() const noexcept { return region_; }
    GenomicRegion& region() noexcept { return region_; }

private:
    GenomicRegion region_;
    NucleotideSequence motif_;
};

unsigned count_periods(const TandemRepeat& repeat) noexcept;

bool operator==(const TandemRepeat& lhs, const TandemRepeat& rhs) noexcept;
bool operator<(const TandemRepeat& lhs, const TandemRepeat& rhs) noexcept;

} // namespace octopus

#endif
