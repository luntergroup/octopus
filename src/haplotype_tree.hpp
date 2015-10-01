//
//  haplotype_tree.hpp
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
#include <boost/graph/adjacency_list.hpp>

#include "allele.hpp"
#include "variant.hpp"
#include "haplotype.hpp"

class GenomicRegion;
class ReferenceGenome;

namespace Octopus
{

class HaplotypeTree
{
public:
    HaplotypeTree() = delete;
    explicit HaplotypeTree(ReferenceGenome& reference);
    ~HaplotypeTree() = default;
    
    HaplotypeTree(const HaplotypeTree&)            = default;
    HaplotypeTree& operator=(const HaplotypeTree&) = default;
    HaplotypeTree(HaplotypeTree&&)                 = default;
    HaplotypeTree& operator=(HaplotypeTree&&)      = default;
    
    bool empty() const noexcept;
    unsigned num_haplotypes() const noexcept;
    bool contains(const Haplotype& haplotype) const;
    bool is_unique(const Haplotype& haplotype) const;
    GenomicRegion get_seperation_region(const Haplotype& first, const Haplotype& second) const;
    void extend(const Allele& allele);
    std::vector<Haplotype> get_haplotypes(const GenomicRegion& region);
    void prune_all(const Haplotype& haplotype);
    void prune_unique(const Haplotype& haplotype);
    void clear();
    
private:
    using Tree   = boost::adjacency_list< boost::listS, boost::listS, boost::bidirectionalS, Allele, boost::no_property>;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree tree_;
    Vertex root_;
    std::list<Vertex> haplotype_leafs_;
    ReferenceGenome& reference_;
    std::unordered_multimap<Haplotype, Vertex> haplotype_leaf_cache_;
    std::unordered_set<Haplotype> recently_removed_haplotypes_;
    
    using LeafIterator  = decltype(haplotype_leafs_)::const_iterator;
    using CacheIterator = decltype(haplotype_leaf_cache_)::iterator;
    
    Vertex get_previous_allele(Vertex allele) const;
    bool allele_exists(Vertex leaf, const Allele& allele) const;
    LeafIterator extend_haplotype(LeafIterator haplotype_leaf, const Allele& new_allele);
    Haplotype get_haplotype(Vertex haplotype_end, const GenomicRegion& region) const;
    
    bool is_branch_exact_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const;
    bool is_branch_equal_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const;
    LeafIterator find_exact_haplotype_leaf(LeafIterator first, LeafIterator last, const Haplotype& haplotype) const;
    LeafIterator find_equal_haplotype_leaf(LeafIterator first, LeafIterator last, const Haplotype& haplotype) const;
    std::pair<Vertex, bool> prune_branch(Vertex leaf, const GenomicRegion& region);
};

// non-member methods

namespace detail
{
    template <typename Container>
    void extend_tree(const Container& container, HaplotypeTree& tree, Allele)
    {
        for (const auto& allele : container) {
            tree.extend(allele);
        }
    }
    
    template <typename Container>
    void extend_tree(const Container& container, HaplotypeTree& tree, Variant)
    {
        for (const auto& variant : container) {
            tree.extend(variant.get_reference_allele());
            tree.extend(variant.get_alternative_allele());
        }
    }
    
} // namespace detail

template <typename Container>
void extend_tree(const Container& container, HaplotypeTree& tree)
{
    detail::extend_tree(container, tree, typename Container::value_type {});
}
    
} // namespace Octopus

#endif /* defined(__Octopus__haplotype_tree__) */
