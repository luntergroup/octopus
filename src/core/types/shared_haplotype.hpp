// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef shared_haplotype_hpp
#define shared_haplotype_hpp

#include <memory>
#include <utility>
#include <type_traits>
#include <functional>

#include "concepts/mappable.hpp"
#include "concepts/comparable.hpp"
#include "haplotype.hpp"

namespace octopus {

class SharedHaplotype : public Mappable<SharedHaplotype> , public Comparable<SharedHaplotype>
{
public:
    using MappingDomain = Haplotype::MappingDomain;
    
    SharedHaplotype() = delete;
    
    SharedHaplotype(const Haplotype& haplotype)
    : haplotype_ {std::make_shared<Haplotype>(haplotype)}
    {}
    SharedHaplotype(Haplotype&& haplotype)
    : haplotype_ {std::make_shared<Haplotype>(std::move(haplotype))}
    {}
    SharedHaplotype(const std::shared_ptr<Haplotype>& haplotype)
    : haplotype_ {haplotype}
    {}
    
    SharedHaplotype(const SharedHaplotype&)            = default;
    SharedHaplotype& operator=(const SharedHaplotype&) = default;
    SharedHaplotype(SharedHaplotype&&)                 = default;
    SharedHaplotype& operator=(SharedHaplotype&&)      = default;
    
    ~SharedHaplotype() = default;
    
    Haplotype& haplotype() noexcept { return *haplotype_; }
    const Haplotype& haplotype() const noexcept { return *haplotype_; }
    
    operator Haplotype&() noexcept { return haplotype(); }
    operator const Haplotype&() const noexcept { return haplotype(); }
    
    decltype(auto) mapped_region() const noexcept { return haplotype().mapped_region(); }
    
    decltype(auto) contains(const ContigAllele& allele) const { return haplotype().contains(allele); }
    decltype(auto) contains(const Allele& allele) const  { return haplotype().contains(allele); }
    
    decltype(auto) includes(const ContigAllele& allele) const { return haplotype().includes(allele); }
    decltype(auto) includes(const Allele& allele) const  { return haplotype().includes(allele); }
    
    decltype(auto) sequence(const ContigRegion& region) const { return haplotype().sequence(region); }
    decltype(auto) sequence(const GenomicRegion& region) const { return haplotype().sequence(region); }
    decltype(auto) sequence() const noexcept { return haplotype().sequence(); }
    
    decltype(auto) sequence_size(const ContigRegion& region) const { return haplotype().sequence_size(region); }
    decltype(auto) sequence_size(const GenomicRegion& region) const { return haplotype().sequence_size(region); }
    
    decltype(auto) difference(const SharedHaplotype& other) const { return haplotype().difference(other.haplotype()); }
    decltype(auto) cigar() const { return haplotype().cigar(); }
    
    friend bool operator==(const SharedHaplotype& lhs, const SharedHaplotype& rhs) noexcept
    {
        return lhs.haplotype_ == rhs.haplotype_ || lhs.haplotype() == rhs.haplotype();
    }
    
private:
    std::shared_ptr<Haplotype> haplotype_;
};

inline bool operator<(const SharedHaplotype& lhs, const SharedHaplotype& rhs)
{
    return lhs.haplotype() < rhs.haplotype();
}

template <typename MappableType, typename HaplotypeType>
std::enable_if_t<std::is_same<MappableType, SharedHaplotype>::value, SharedHaplotype>
copy(const HaplotypeType& haplotype, const GenomicRegion& region)
{
    return {std::make_shared<Haplotype>(copy<Haplotype>(haplotype.haplotype(), region))};
}

} // namespace octopus

namespace std {

template <>
struct hash<octopus::SharedHaplotype>
{
    size_t operator()(const octopus::SharedHaplotype& haplotpe) const
    {
        return hash<octopus::Haplotype>{}(haplotpe.haplotype());
    }
};

} // namespace std

namespace boost {

template <>
struct hash<octopus::SharedHaplotype>
{
    std::size_t operator()(const octopus::SharedHaplotype& h) const
    {
        return std::hash<octopus::SharedHaplotype>{}(h);
    }
};

} // namespace boost

#endif
