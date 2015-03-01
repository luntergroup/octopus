//
//  variant_candidate_generator.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_candidate_generator__
#define __Octopus__variant_candidate_generator__

#include <set>
#include <vector>
#include <memory> // std::unique_ptr
#include <algorithm> // std::for_each
#include <iterator> // std::make_move_iterator

#include "i_variant_candidate_generator.h"
#include "reference_genome.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "variant.h"

//TODO: profile to see if std::set is useful here. Note we cannot std::move from
// a std::set so we cannot 'move' merge results.. maybe std::vector would be better.

class VariantCandidateGenerator : public IVariantCandidateGenerator
{
public:
    VariantCandidateGenerator() = default;
    ~VariantCandidateGenerator() override = default;
    
    VariantCandidateGenerator(const VariantCandidateGenerator&)            = default;
    VariantCandidateGenerator& operator=(const VariantCandidateGenerator&) = default;
    VariantCandidateGenerator(VariantCandidateGenerator&&)                 = default;
    VariantCandidateGenerator& operator=(VariantCandidateGenerator&&)      = default;
    
    void register_generator(std::unique_ptr<IVariantCandidateGenerator> g);
    void add_read(const AlignedRead& a_read) override;
    std::set<Variant> get_candidates(const GenomicRegion& a_region) override;
    void clear() override;
    template <typename ReadIterator> void add_reads(ReadIterator begin, ReadIterator end);
    
private:
    std::vector<std::unique_ptr<IVariantCandidateGenerator>> generator_list_;
};

inline void VariantCandidateGenerator::register_generator(std::unique_ptr<IVariantCandidateGenerator> g)
{
    generator_list_.emplace_back(std::move(g));
}

inline void VariantCandidateGenerator::add_read(const AlignedRead& a_read)
{
    std::for_each(std::begin(generator_list_), std::end(generator_list_), [&a_read] (auto& g) {
        g->add_read(a_read);
    });
}

inline std::set<Variant> VariantCandidateGenerator::get_candidates(const GenomicRegion& a_region)
{
    std::set<Variant> result {};
    std::for_each(std::begin(generator_list_), std::end(generator_list_), [&result, &a_region] (auto& g) {
        auto generator_result = g->get_candidates(a_region);
        result.insert(std::begin(generator_result),
                      std::end(generator_result));
    });
    return result;
}

inline void VariantCandidateGenerator::clear()
{
    std::for_each(std::begin(generator_list_), std::end(generator_list_), [] (auto& g) {
        g->clear();
    });
}

template <typename ReadIterator>
void VariantCandidateGenerator::add_reads(ReadIterator begin, ReadIterator end)
{
    std::for_each(begin, end, [this] (const auto& a_read ) { add_read(a_read); });
}

#endif /* defined(__Octopus__variant_candidate_generator__) */
