//
//  variant.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_h
#define Octopus_variant_h

#include <string>
#include <cstdint>
#include <cstddef>
#include <functional>

#include "genomic_region.h"
#include "comparable.h"

using std::uint_fast32_t;
using std::size_t;

/**
   All variants are considered replacements, i.e. a region of the reference contig is removed,
   and replaced with a novel sequence.
 
   The underlying logic is designed to be transparent to the variant type. The exception is
   modelling the prior probability of different variants. This solution uses a strategy pattern
   injection (the functional 'prior_model' is injected) to acheive runtime polymorphism.
 */
class Variant : Comparable<Variant>
{
public:
    using SizeType = GenomicRegion::SizeType;
    
    Variant() = delete;
    template <typename T1, typename T2, typename T3>
    explicit Variant(T1&& reference_contig_name, SizeType reference_removed_region_begin,
                     T2&& reference_sequence_removed, T3&& sequence_added,
                     std::function<double()> prior_model);
    template <typename T1, typename T2>
    explicit Variant(GenomicRegion reference_removed_region, T1&& reference_sequence_removed,
                     T2&& sequence_added, std::function<double()> prior_model);
    
    Variant(const Variant&)            = default;
    Variant& operator=(const Variant&) = default;
    Variant(Variant&&)                 = default;
    Variant& operator=(Variant&&)      = default;
    
    const GenomicRegion& get_removed_region() const noexcept;
    SizeType get_removed_region_begin() const noexcept;
    SizeType get_removed_region_end() const noexcept;
    const std::string& get_sequence_added() const noexcept;
    const std::string& get_sequence_removed() const noexcept;
    bool has_support(unsigned long num_reads) const noexcept;
    unsigned long get_num_supporting_reads() const noexcept;
    double get_prior_probability() const noexcept;
    void add_support(unsigned long num_reads) noexcept;
    
private:
    GenomicRegion reference_removed_region_;
    std::string reference_sequence_removed_;
    std::string sequence_added_;
    unsigned long num_supporting_reads_;
    std::function<double()> prior_model_;
};

template <typename T1, typename T2, typename T3>
Variant::Variant(T1&& reference_contig_name, SizeType reference_removed_region_begin,
                 T2&& reference_sequence_removed, T3&& sequence_added,
                 std::function<double()> prior_model)
:reference_removed_region_ {reference_contig_name, reference_removed_region_begin,
                           reference_removed_region_begin +
    static_cast<SizeType>(reference_sequence_removed.size())},
reference_sequence_removed_ {std::forward<T2>(reference_sequence_removed)},
sequence_added_ {std::forward<T3>(sequence_added)},
prior_model_ {prior_model}
{}

template <typename T1, typename T2>
Variant::Variant(GenomicRegion reference_removed_region, T1&& reference_sequence_removed,
                 T2&& sequence_added, std::function<double()> prior_model)
:reference_removed_region_ {std::move(reference_removed_region)},
reference_sequence_removed_ {std::forward<T1>(reference_sequence_removed)},
sequence_added_ {std::forward<T2>(sequence_added)},
prior_model_ {prior_model}
{}

inline const GenomicRegion& Variant::get_removed_region() const noexcept
{
    return reference_removed_region_;
}

inline Variant::SizeType Variant::get_removed_region_begin() const noexcept
{
    return reference_removed_region_.get_begin();
}

inline Variant::SizeType Variant::get_removed_region_end() const noexcept
{
    return reference_removed_region_.get_end();
}

inline const std::string& Variant::get_sequence_added() const noexcept
{
    return sequence_added_;
}

inline const std::string& Variant::get_sequence_removed() const noexcept
{
    return reference_sequence_removed_;
}

inline bool Variant::has_support(unsigned long num_reads) const noexcept
{
    return num_supporting_reads_ >= num_reads;
}

inline unsigned long Variant::get_num_supporting_reads() const noexcept
{
    return num_supporting_reads_;
}

inline void Variant::add_support(unsigned long num_reads) noexcept
{
    num_supporting_reads_ += num_reads;
}

inline double Variant::get_prior_probability() const noexcept
{
    return prior_model_();
}

inline bool overlaps(const Variant& lhs, const Variant& rhs) noexcept
{
    return overlaps(lhs.get_removed_region(), rhs.get_removed_region());
}

inline bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_removed_region() == rhs.get_removed_region() &&
           lhs.get_sequence_added() == rhs.get_sequence_added() &&
           lhs.get_sequence_removed() == rhs.get_sequence_removed();
}

inline bool operator<(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_removed_region() < rhs.get_removed_region();
}

namespace std {
    template <> struct hash<Variant>
    {
        size_t operator()(const Variant& v) const
        {
            return hash<GenomicRegion>()(v.get_removed_region());
        }
    };
}

#endif
