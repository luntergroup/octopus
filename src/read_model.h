//
//  read_model.h
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_model__
#define __Octopus__read_model__

#include <vector>
#include <string>
#include <unordered_map>
#include <cstddef> // std::size_t
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "common.h"
#include "haplotype.h"
#include "genotype.h"
#include "aligned_read.h"

namespace Octopus
{

class ReadModel
{
public:
    using RealType     = Octopus::ProbabilityType;
    using SampleIdType = Octopus::SampleIdType;
    
    ReadModel() = delete;
    explicit ReadModel(unsigned ploidy, bool can_cache_reads = true);
    ~ReadModel() = default;
    
    ReadModel(const ReadModel&)            = default;
    ReadModel& operator=(const ReadModel&) = default;
    ReadModel(ReadModel&&)                 = default;
    ReadModel& operator=(ReadModel&&)      = default;
    
    // ln p(read | haplotype)
    RealType log_probability(const AlignedRead& read, const Haplotype& haplotype, SampleIdType sample);
    
    // ln p(read | genotype)
    RealType log_probability(const AlignedRead& read, const Genotype<Haplotype>& genotype, SampleIdType sample);
    
    // ln p(reads | genotype)
    template <typename ForwardIterator>
    RealType log_probability(ForwardIterator first_read, ForwardIterator last_read,
                             const Genotype<Haplotype>& genotype, SampleIdType sample);
    
    void clear_cache();
    
private:
    unsigned ploidy_;
    bool can_cache_reads_;
    std::unordered_map<SampleIdType, std::unordered_map<AlignedRead,
                        std::unordered_map<Haplotype, RealType>>> read_log_probability_cache_;
    
    using GenotypeReadKey = std::tuple<Genotype<Haplotype>, AlignedRead, std::size_t>;
    
    struct GenotypeReadKeyHash
    {
        std::size_t operator()(const GenotypeReadKey& key) const;
    };
    
    std::unordered_map<SampleIdType, std::unordered_map<GenotypeReadKey, RealType, GenotypeReadKeyHash>>
        genotype_log_probability_cache_;
    
    const double ln_ploidy_;
    
    bool is_read_in_cache(SampleIdType sample, const AlignedRead& read, const Haplotype& haplotype) const noexcept;
    
    void add_read_to_cache(SampleIdType sample, const AlignedRead& read, const Haplotype& haplotype,
                           RealType read_log_probability);
    
    RealType get_log_probability_from_cache(SampleIdType sample, const AlignedRead& read,
                                            const Haplotype& haplotype) const;
    
    template <typename ForwardIterator>
    bool is_genotype_in_cache(SampleIdType sample, ForwardIterator first_read, ForwardIterator last_read,
                              const Genotype<Haplotype>& genotype) const noexcept;
    
    template <typename ForwardIterator>
    void add_genotype_to_cache(SampleIdType sample, ForwardIterator first_read, ForwardIterator last_read,
                               const Genotype<Haplotype>& genotype, RealType genotype_log_probability);
    
    template <typename ForwardIterator>
    RealType get_log_probability_from_cache(SampleIdType sample, ForwardIterator first_read, ForwardIterator last_read,
                                            const Genotype<Haplotype>& genotype) const;
    
    // These are just for optimisation
    RealType log_probability_haploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                     SampleIdType sample);
    RealType log_probability_diploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                     SampleIdType sample);
    RealType log_probability_triploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                      SampleIdType sample);
    RealType log_probability_polyploid(const AlignedRead& read, const Genotype<Haplotype>& genotype,
                                       SampleIdType sample);
};

// ln p(reads | genotype) = sum (read in reads} ln p(read | genotype)
template <typename ForwardIterator>
ReadModel::RealType ReadModel::log_probability(ForwardIterator first_read, ForwardIterator last_read,
                                               const Genotype<Haplotype>& genotype, SampleIdType sample)
{
    if (is_genotype_in_cache(sample, first_read, last_read, genotype)) {
        return get_log_probability_from_cache(sample, first_read, last_read, genotype);
    }
    
    RealType result {0};
    
    std::for_each(first_read, last_read, [this, &genotype, &sample, &result] (const auto& read) {
        result += log_probability(read, genotype, sample);
    });
    
    add_genotype_to_cache(sample, first_read, last_read, genotype, result);
    
    return result;
}

template <typename ForwardIterator>
inline
bool ReadModel::is_genotype_in_cache(SampleIdType sample, ForwardIterator first_read, ForwardIterator last_read,
                                     const Genotype<Haplotype>& genotype) const noexcept
{
    if (genotype_log_probability_cache_.count(sample) == 0) return false;
    return genotype_log_probability_cache_.at(sample).count(std::make_tuple(genotype, *first_read,
                                                                            std::distance(first_read, last_read))) > 0;
}

template <typename ForwardIterator>
inline
void ReadModel::add_genotype_to_cache(SampleIdType sample, ForwardIterator first_read, ForwardIterator last_read,
                                      const Genotype<Haplotype>& genotype, RealType genotype_log_probability)
{
    genotype_log_probability_cache_[sample][std::make_tuple(genotype, *first_read,
                                                            std::distance(first_read, last_read))] = genotype_log_probability;
}

template <typename ForwardIterator>
inline
ReadModel::RealType
ReadModel::get_log_probability_from_cache(SampleIdType sample, ForwardIterator first_read, ForwardIterator last_read,
                                          const Genotype<Haplotype>& genotype) const
{
    return genotype_log_probability_cache_.at(sample).at(std::make_tuple(genotype, *first_read,
                                                                         std::distance(first_read, last_read)));
}

inline std::size_t ReadModel::GenotypeReadKeyHash::operator()(const GenotypeReadKey &key) const
{
    std::size_t seed {};
    boost::hash_combine(seed, std::hash<Genotype<Haplotype>>()(std::get<0>(key)));
    boost::hash_combine(seed, std::hash<AlignedRead>()(std::get<1>(key)));
    boost::hash_combine(seed, std::get<2>(key));
    return seed;
}

} // end namespace Octopus

#endif /* defined(__Octopus__read_model__) */
