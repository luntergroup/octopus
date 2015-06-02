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
#include <cstddef> // std::size_t
#include <boost/graph/adjacency_list.hpp>

#include "allele.h"
#include "haplotype.h"

class GenomicRegion;
class ReferenceGenome;

namespace Octopus
{

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
    
    // Non-modifying
    bool empty() const noexcept;
    std::size_t num_haplotypes() const noexcept;
    bool contains(const Haplotype& a_haplotype) const;
    bool is_unique(const Haplotype& a_haplotype) const;
    GenomicRegion get_seperation_region(const Haplotype& first, const Haplotype& second) const;
    void extend(const Allele& an_allele);
    Haplotypes get_haplotypes(const GenomicRegion& a_region);
    
    // Modifying
    void prune_all(const Haplotype& a_haplotype);
    void prune_unique(const Haplotype& a_haplotype);
    void clear();
    
private:
    using Tree = boost::adjacency_list<
        boost::listS, boost::listS, boost::bidirectionalS, Allele, boost::no_property
    >;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree the_tree_;
    Vertex the_root_;
    std::list<Vertex> the_haplotype_leafs_;
    ReferenceGenome& the_reference_;
    std::unordered_multimap<Haplotype, Vertex> haplotype_leaf_cache_;
    std::unordered_set<Haplotype> recently_removed_haplotypes_;
    
    using LeafIterator = decltype(the_haplotype_leafs_)::const_iterator;
    using CacheIterator = decltype(haplotype_leaf_cache_)::iterator;
    
    Vertex get_previous_allele(Vertex an_allele) const;
    bool allele_exists(Vertex allele, const Allele& an_allele) const;
    LeafIterator extend_haplotype(LeafIterator haplotype_leaf, const Allele& the_new_allele);
    Haplotype get_haplotype(Vertex haplotype_end, const GenomicRegion& a_region) const;
    
    bool is_branch_exact_haplotype(Vertex haplotype_end, const Haplotype& a_haplotype) const;
    bool is_branch_equal_haplotype(Vertex haplotype_end, const Haplotype& a_haplotype) const;
    LeafIterator find_exact_haplotype_leaf(LeafIterator first, LeafIterator last, const Haplotype& a_haplotype) const;
    LeafIterator find_equal_haplotype_leaf(LeafIterator first, LeafIterator last, const Haplotype& a_haplotype) const;
    std::pair<Vertex, bool> prune_branch(Vertex leaf, const GenomicRegion& a_region);
};

} // end namespace Octopus

#endif /* defined(__Octopus__haplotype_tree__) */
