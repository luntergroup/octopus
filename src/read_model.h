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

#include "haplotype.h"
#include "genotype.h"
#include "aligned_read.h"

class ReadModel
{
public:
    using RealType     = double;
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    using SampleIdType = std::string;
    
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
    RealType log_probability(const AlignedRead& read, const Genotype& genotype, SampleIdType sample);
    
    // ln p(reads | genotype)
    RealType log_probability(ReadIterator first, ReadIterator last, const Genotype& genotype, SampleIdType sample);
    
    void clear_cache();
    
private:
    unsigned ploidy_;
    bool can_cache_reads_;
    std::unordered_map<SampleIdType, std::unordered_map<AlignedRead,
                        std::unordered_map<Haplotype, RealType>>> read_log_probability_cache_;
    
    using GenotypeReadKey = std::tuple<Genotype, AlignedRead, std::size_t>;
    
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
    
    bool is_genotype_in_cache(SampleIdType sample, ReadIterator first, ReadIterator last,
                              const Genotype& genotype) const noexcept;
    
    void add_genotype_to_cache(SampleIdType sample, ReadIterator first, ReadIterator last,
                               const Genotype& genotype, RealType genotype_log_probability);
    
    RealType get_log_probability_from_cache(SampleIdType sample, ReadIterator first, ReadIterator last,
                                            const Genotype& genotype) const;
    
    // These are just for optimisation
    RealType log_probability_haploid(const AlignedRead& read, const Genotype& genotype,
                                     SampleIdType sample);
    RealType log_probability_diploid(const AlignedRead& read, const Genotype& genotype,
                                     SampleIdType sample);
    RealType log_probability_triploid(const AlignedRead& read, const Genotype& genotype,
                                      SampleIdType sample);
    RealType log_probability_polyploid(const AlignedRead& read, const Genotype& genotype,
                                       SampleIdType sample);
};

inline std::size_t ReadModel::GenotypeReadKeyHash::operator()(const GenotypeReadKey &key) const
{
    std::size_t seed {};
    boost::hash_combine(seed, std::hash<Genotype>()(std::get<0>(key)));
    boost::hash_combine(seed, std::hash<AlignedRead>()(std::get<1>(key)));
    boost::hash_combine(seed, std::get<2>(key));
    return seed;
}

#endif /* defined(__Octopus__read_model__) */
