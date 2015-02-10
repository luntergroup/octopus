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

using std::uint_fast32_t;
using std::size_t;

/**
   All variants are considered replacements, e.g. a SNP can be represented
   as a replacement of one nucleotide with another.
 
   The underlying logic is designed to be transparent to the variant type. The exception is
   modelling the prior probability of different variants. This solution uses a strategy pattern
   injection (the functional 'prior_model' is injected) to acheive runtime polymorphism.
 */
class Variant
{
public:
    
    Variant() = delete;
    Variant(std::string contig_name, int_fast32_t contig_start_pos, std::string sequence_added,
            std::string sequence_removed, std::function<double()> prior_model);
    ~Variant();
    
    GenomicRegion get_region() const noexcept;
    const std::string& get_sequence_added() const noexcept;
    const std::string& get_sequence_removed() const noexcept;
    unsigned long get_num_supporting_reads() const noexcept;
    double get_prior_probability() const noexcept;
    
    void add_support(unsigned long num_reads) noexcept;
    
private:
    GenomicRegion contig_region_;
    std::string sequence_added_;
    std::string sequence_removed_;
    unsigned long num_supporting_reads_;
    std::function<double()> prior_model_;
};

inline
Variant::Variant(std::string contig_name, int_fast32_t contig_start_pos, std::string sequence_added,
                 std::string sequence_removed, std::function<double()> prior_model)
: contig_region_(contig_name, contig_start_pos, contig_start_pos + sequence_removed.size()),
  sequence_added_(sequence_added),
  sequence_removed_(sequence_removed),
  prior_model_(prior_model)
{}

inline Variant::~Variant() {}

inline GenomicRegion Variant::get_region() const noexcept
{
    return contig_region_;
}

inline const std::string& Variant::get_sequence_added() const noexcept
{
    return sequence_added_;
}

inline const std::string& Variant::get_sequence_removed() const noexcept
{
    return sequence_removed_;
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
    return overlaps(lhs.get_region(), rhs.get_region());
}

inline bool operator==(const Variant& lhs, const Variant& rhs)
{
    return lhs.get_region() == rhs.get_region() &&
    lhs.get_sequence_added() == rhs.get_sequence_added() &&
    lhs.get_sequence_removed() == rhs.get_sequence_removed();
}
inline bool operator!=(const Variant& lhs, const Variant& rhs) {return !operator==(lhs,rhs);}

namespace std {
    template <> struct hash<Variant>
    {
        size_t operator()(const Variant& v) const
        {
            return hash<size_t>()(v.get_region().get_begin_pos()); //TODO: do something better!
        }
    };
}

#endif
