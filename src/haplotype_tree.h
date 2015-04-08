//
//  haplotype_tree.h
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype_tree__
#define __Octopus__haplotype_tree__

#include <vector>
#include <list>
#include <cstddef>
#include <boost/graph/adjacency_list.hpp>

#include "haplotype.h"
#include "read_manager.h"
#include "genomic_region.h"

class ReferenceGenome;
class VariationalBayesGenotypeModel;
class Variant;

using std::size_t;

class HaplotypeTree
{
public:
    using Haplotypes = std::vector<Haplotype>;
    
    HaplotypeTree() = delete;
    explicit HaplotypeTree(ReferenceGenome& the_reference, ReadManager& the_reads,
                           VariationalBayesGenotypeModel& the_genotype_model,
                           const std::vector<ReadManager::SampleIdType>& the_sample_ids,
                           size_t max_num_haplotypes, double min_posterior);
    ~HaplotypeTree() = default;
    
    HaplotypeTree(const HaplotypeTree&)            = default;
    HaplotypeTree& operator=(const HaplotypeTree&) = default;
    HaplotypeTree(HaplotypeTree&&)                 = default;
    HaplotypeTree& operator=(HaplotypeTree&&)      = default;
    
    Haplotypes get_haplotypes(const std::vector<Variant>& ordered_variants);
    
private:
    struct Allele
    {
        Allele() = default;
        template <typename T>
        Allele(const GenomicRegion& the_reference_region, T&& the_sequence)
        :
        the_reference_region {the_reference_region},
        the_sequence {std::forward<T>(the_sequence)}
        {}
        
        GenomicRegion the_reference_region;
        Variant::SequenceType the_sequence;
    };
    
    struct HaplotypeNode
    {
        Allele the_variant;
        double haplotype_probability;
    };
    
    using Tree = boost::adjacency_list<
        boost::listS, boost::listS, boost::undirectedS, HaplotypeNode, boost::no_property
    >;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    ReferenceGenome& the_reference_;
    ReadManager& the_reads_;
    VariationalBayesGenotypeModel& the_genotype_model_;
    std::vector<ReadManager::SampleIdType> the_sample_ids_;
    
    size_t max_num_haplotypes_;
    double min_posterior_;
    size_t num_extensions_since_posterior_update_;
    
    Tree the_tree_;
    std::list<Vertex> haplotype_branch_ends_;
    
    void init_tree();
    size_t num_haplotypes() const;
    void extend_tree(const Variant& a_variant);
    std::vector<Vertex> extend_haplotype(Vertex haplotype_branch_end, const Variant& a_variant);
    Haplotypes get_unupdated_haplotypes();
    void update_posteriors();
    void prune_low_probability_haplotypes(size_t n);
    void prune_haplotype(Vertex haplotype_end);
    Haplotypes get_highest_probability_haplotypes();
};

#endif /* defined(__Octopus__haplotype_tree__) */
