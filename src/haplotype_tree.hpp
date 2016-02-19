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
#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <type_traits>

#include <boost/graph/adjacency_list.hpp>

#include "variant.hpp"
#include "haplotype.hpp"
#include "genomic_region.hpp"
#include "reference_genome.hpp"

namespace Octopus
{

class HaplotypeTree
{
public:
    using ContigNameType = GenomicRegion::ContigNameType;
    
    HaplotypeTree() = delete;
    explicit HaplotypeTree(const ContigNameType& contig, const ReferenceGenome& reference);
    ~HaplotypeTree() = default;
    
    HaplotypeTree(const HaplotypeTree&)            = default;
    HaplotypeTree& operator=(const HaplotypeTree&) = default;
    HaplotypeTree(HaplotypeTree&&)                 = default;
    HaplotypeTree& operator=(HaplotypeTree&&)      = default;
    
    bool empty() const noexcept;
    unsigned num_haplotypes() const noexcept;
    bool contains(const Haplotype& haplotype) const;
    bool is_unique(const Haplotype& haplotype) const;
    HaplotypeTree& extend(const ContigAllele& allele);
    HaplotypeTree& extend(const Allele& allele);
    
    GenomicRegion get_region() const;
    
    std::vector<Haplotype> extract_haplotypes() const;
    std::vector<Haplotype> extract_haplotypes(const GenomicRegion& region) const;
    
    void prune_all(const Haplotype& haplotype);
    void prune_unique(const Haplotype& haplotype);
    void remove(const GenomicRegion& region);
    void clear();
    
private:
    using Tree = boost::adjacency_list<
                     boost::listS, boost::listS, boost::bidirectionalS,
                     ContigAllele, boost::no_property
                 >;
    
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    
    Tree tree_;
    Vertex root_;
    std::list<Vertex> haplotype_leafs_;
    
    GenomicRegion region_;
    
    mutable std::unordered_multimap<Haplotype, Vertex> haplotype_leaf_cache_;
    mutable std::unordered_set<Haplotype> recently_removed_haplotypes_;
    
    using LeafIterator  = decltype(haplotype_leafs_)::const_iterator;
    using CacheIterator = decltype(haplotype_leaf_cache_)::iterator;
    
    bool is_bifurcating(Vertex v) const;
    Vertex remove_forward(Vertex u);
    Vertex remove_backward(Vertex v);
    Vertex get_previous_allele(Vertex allele) const;
    bool allele_exists(Vertex leaf, const ContigAllele& allele) const;
    LeafIterator extend_haplotype(LeafIterator leaf, const ContigAllele& new_allele);
    Haplotype extract_haplotype(Vertex leaf, const GenomicRegion& region) const;
    bool define_same_haplotype(Vertex leaf1, Vertex leaf2) const;
    bool is_branch_exact_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const;
    bool is_branch_equal_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const;
    LeafIterator find_exact_haplotype_leaf(LeafIterator first, LeafIterator last,
                                           const Haplotype& haplotype) const;
    LeafIterator find_equal_haplotype_leaf(LeafIterator first, LeafIterator last,
                                           const Haplotype& haplotype) const;
    std::pair<Vertex, bool> remove(Vertex leaf, const ContigRegion& region);
    std::pair<Vertex, bool> remove_external(Vertex leaf, const ContigRegion& region);
    std::pair<Vertex, bool> remove_internal(Vertex leaf, const ContigRegion& region);
};

// non-member methods

namespace detail
{
    template <typename InputIt>
    void extend_tree(InputIt first, InputIt last, HaplotypeTree& tree, Allele)
    {
        std::for_each(first, last, [&] (const auto& allele) {
            tree.extend(allele);
        });
    }
    
    template <typename InputIt>
    void extend_tree(InputIt first, InputIt last, HaplotypeTree& tree, Variant)
    {
        std::for_each(first, last, [&] (const auto& variant) {
            tree.extend(variant.get_ref_allele());
            tree.extend(variant.get_alt_allele());
        });
    }
    
    template <typename T>
    constexpr bool is_variant_or_allele = std::is_same<T, Allele>::value || std::is_same<T, Variant>::value;
} // namespace detail

template <typename InputIt>
void extend_tree(InputIt first, InputIt last, HaplotypeTree& tree)
{
    using MappableType = std::decay_t<typename std::iterator_traits<InputIt>::value_type>;
    
    static_assert(detail::is_variant_or_allele<MappableType>,
                  "extend_tree only works for containers of Alleles or Variants");
    
    detail::extend_tree(first, last, tree, MappableType {});
}

template <typename Container>
void extend_tree(const Container& elements, HaplotypeTree& tree)
{
    extend_tree(std::cbegin(elements), std::cend(elements), tree);
}

template <typename InputIt>
std::vector<Haplotype>
generate_all_haplotypes(InputIt first, InputIt last, const ReferenceGenome& reference)
{
    using MappableType = std::decay_t<typename std::iterator_traits<InputIt>::value_type>;
    
    static_assert(detail::is_variant_or_allele<MappableType>,
                  "generate_all_haplotypes only works for containers of Alleles or Variants");
    
    if (first == last) return {};
    
    HaplotypeTree tree {contig_name(*first), reference};
    
    extend_tree(first, last, tree);
    
    return tree.extract_haplotypes();
}

template <typename Container>
auto generate_all_haplotypes(const Container& elements, const ReferenceGenome& reference)
{
    return generate_all_haplotypes(std::cbegin(elements), std::cend(elements), reference);
}
    
} // namespace Octopus

#endif /* defined(__Octopus__haplotype_tree__) */
