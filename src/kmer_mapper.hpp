//
//  kmer_mapper.hpp
//  Octopus
//
//  Created by Daniel Cooke on 23/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef kmer_mapper_hpp
#define kmer_mapper_hpp

#include <vector>
#include <unordered_map>
#include <cstddef>
#include <functional>

#include "kmer_mapping.hpp"
#include "common.hpp"
#include "haplotype.hpp"
#include "aligned_read.hpp"

namespace Octopus
{
    class KmerMapper
    {
    public:
        KmerMapper() = delete;
        explicit KmerMapper(const ReadMap& reads, const std::vector<Haplotype>& haplotypes);
        ~KmerMapper() = default;
        
        std::vector<std::size_t> map(const AlignedRead& read, const Haplotype& haplotype) const;
        
        void clear() noexcept;
        
    private:
        static constexpr unsigned char KMER_SIZE {5};
        
        using ReadKeyType      = AlignedRead::SequenceType;
        using HaplotypeKeyType = Haplotype::SequenceType;
        
        mutable std::unordered_map<ReadKeyType, KmerPerfectHashes> read_cache_;
        mutable std::unordered_map<ReadKeyType, KmerHashTable> haplotype_cache_;
    };
}

#endif /* kmer_mapper_hpp */
