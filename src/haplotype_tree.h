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
#include <unordered_map>
#include <unordered_set>
#include <cstddef>
#include <boost/graph/adjacency_list.hpp>

#include "allele.h"
#include "haplotype.h"
#include "genotype.h"
#include "genomic_region.h"

class ReferenceGenome;
class Variant;

class HaplotypeTree
{
public:
    using Haplotypes = std::vector<Haplotype>;
    
    HaplotypeTree() = delete;
    explicit HaplotypeTree(ReferenceGenome& the_reference);
    ~HaplotypeTree() = default;
    
    HaplotypeTree(const HaplotypeTree&)            = default;
    HaplotypeTree& operator=(const HaplotypeTree&) = default;
    HaplotypeTree(HaplotypeTree&&)                 = default;
    HaplotypeTree& operator=(HaplotypeTree&&)      = default;
    
    std::size_t num_haplotypes() const;
    void extend(const Allele& an_allele);
    Haplotypes get_haplotypes(const GenomicRegion& a_region);
    void prune_all(const Haplotype& haplotype);
    void prune_unique(const Haplotype& haplotype);
    
private:
    struct AlleleNode
    {
        Allele the_allele;
    };
    
    using Tree = boost::adjacency_list<
        boost::listS, boost::listS, boost::bidirectionalS, AlleleNode, boost::no_property
    >;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree the_tree_;
    Vertex the_root_;
    std::list<Vertex> haplotype_leafs_;
    ReferenceGenome& the_reference_;
    std::unordered_multimap<Haplotype, Vertex> haplotype_leaf_cache_;
    std::unordered_set<Haplotype> recently_removed_haplotypes_;
    
    using LeafIterator = decltype(haplotype_leafs_)::const_iterator;
    using CacheIterator = decltype(haplotype_leaf_cache_)::iterator;
    
    Vertex get_previous_allele(Vertex allele) const;
    bool allele_exists(Vertex allele, const Allele& an_allele) const;
    LeafIterator extend_haplotype(LeafIterator haplotype_leaf, const Variant& a_variant);
    LeafIterator extend_haplotype(LeafIterator haplotype_leaf, const Allele& the_new_allele);
    Haplotype get_haplotype(Vertex haplotype_end, const GenomicRegion& a_region);
    
    bool is_branch_the_haplotype(Vertex haplotype_end, const Haplotype& haplotype) const;
    LeafIterator find_haplotype_leaf(LeafIterator first, LeafIterator last, const Haplotype& haplotype) const;
    std::pair<Vertex, bool> prune_branch(Vertex leaf, const GenomicRegion& a_region);
    
    void add_to_cache(const Haplotype& haplotype, Vertex leaf);
    std::pair<LeafIterator, bool> get_leaf_from_cache(const Haplotype& haplotype);
};

#endif /* defined(__Octopus__haplotype_tree__) */
