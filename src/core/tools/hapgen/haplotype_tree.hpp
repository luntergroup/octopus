// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef haplotype_tree_hpp
#define haplotype_tree_hpp

#include <vector>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <type_traits>
#include <cmath>
#include <cstddef>

#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>
#include <boost/filesystem.hpp>

#include "basics/genomic_region.hpp"
#include "core/types/allele.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/variant.hpp"
#include "containers/mappable_block.hpp"

namespace octopus {

class ReferenceGenome;

namespace coretools {

class HaplotypeTree
{
public:
    using HaplotypeBlock = MappableBlock<Haplotype>;
    
    using HaplotypeLength = Haplotype::NucleotideSequence::size_type;
    
    HaplotypeTree() = delete;
    
    HaplotypeTree(const GenomicRegion::ContigName& contig, const ReferenceGenome& reference);
    
    HaplotypeTree(const HaplotypeTree&);
    HaplotypeTree& operator=(const HaplotypeTree&);
    HaplotypeTree(HaplotypeTree&&)            = default;
    HaplotypeTree& operator=(HaplotypeTree&&) = default;
    
    ~HaplotypeTree() = default;
    
    bool is_empty() const noexcept;
    
    std::size_t num_haplotypes() const noexcept;
    
    // uses Haplotype::operator== logic
    bool contains(const Haplotype& haplotype) const;
    
    // uses Haplotype::HaveSameAlleles logic
    bool includes(const Haplotype& haplotype) const;
    
    // using Haplotype::HaveSameAlleles logic
    bool is_unique(const Haplotype& haplotype) const;
    
    // Only extends existing leafs
    HaplotypeTree& extend(const ContigAllele& allele);
    HaplotypeTree& extend(const Allele& allele);
    HaplotypeTree& extend(const Haplotype& haplotype);
    
    // Splices into the tree wherever allele can be made a new leaf
    void splice(const ContigAllele& allele);
    void splice(const Allele& allele);
    
    GenomicRegion encompassing_region() const;
    
    HaplotypeBlock extract_haplotypes() const;
    HaplotypeBlock extract_haplotypes(const GenomicRegion& region) const;
    
    std::vector<HaplotypeLength> extract_haplotype_lengths() const;
    std::vector<HaplotypeLength> extract_haplotype_lengths(const GenomicRegion& region) const;
    
    // Using Haplotype::operator== logic
    void prune_all(const Haplotype& haplotype);
    
    // Using Haplotype::HaveSameAlleles logic
    void prune_unique(const Haplotype& haplotype);
    
    void clear(const GenomicRegion& region);
    
    void clear() noexcept;
    
    void write_dot(std::ostream& out) const;
    
private:
    using Tree = boost::adjacency_list<
        boost::listS, boost::listS, boost::bidirectionalS, ContigAllele, boost::no_property
    >;
    
    using Vertex = boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = boost::graph_traits<Tree>::edge_descriptor;
    
    using HaplotypeVertexMultiMap = std::unordered_multimap<Haplotype, Vertex>;
    
    std::reference_wrapper<const ReferenceGenome> reference_;
    Tree tree_;
    Vertex root_;
    std::list<Vertex> haplotype_leafs_;
    GenomicRegion::ContigName contig_;
    
    mutable HaplotypeVertexMultiMap haplotype_leaf_cache_;
    mutable boost::optional<GenomicRegion> tree_region_;
    
    using LeafIterator = decltype(haplotype_leafs_)::iterator;
    using LeafConstIterator = decltype(haplotype_leafs_)::const_iterator;
    using CacheIterator = decltype(haplotype_leaf_cache_)::iterator;
    
    bool is_leaf(Vertex v) const;
    bool is_bifurcating(Vertex v) const;
    Vertex remove_forward(Vertex u);
    Vertex remove_backward(Vertex v);
    Vertex get_previous_allele(Vertex allele) const;
    Vertex find_allele_before(Vertex v, const ContigAllele& allele) const;
    Vertex find_allele_on_branch(LeafIterator leaf, const ContigAllele& allele) const;
    bool allele_exists(Vertex leaf, const ContigAllele& allele) const;
    std::pair<LeafIterator, bool> extend_haplotype(LeafIterator leaf, const ContigAllele& new_allele);
    LeafIterator extend_haplotype(LeafIterator leaf, const Haplotype& other);
    Haplotype extract_haplotype(Vertex leaf, const GenomicRegion& region) const;
    HaplotypeLength extract_haplotype_length(Vertex leaf, const GenomicRegion& region) const;
    bool define_same_haplotype(Vertex leaf1, Vertex leaf2) const;
    bool is_branch_exact_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const;
    bool is_branch_equal_haplotype(Vertex branch_vertex, const Haplotype& haplotype) const;
    LeafConstIterator
     find_exact_haplotype_leaf(LeafConstIterator first, LeafConstIterator last,
                               const Haplotype& haplotype) const;
    LeafConstIterator 
    find_equal_haplotype_leaf(LeafConstIterator first, LeafConstIterator last,
                              const Haplotype& haplotype) const;
    void clear_overlapped(const ContigRegion& region);
    std::pair<Vertex, bool> clear(Vertex leaf, const ContigRegion& region);
    std::pair<Vertex, bool> clear_external(Vertex leaf, const ContigRegion& region);
    std::pair<Vertex, bool> clear_internal(Vertex leaf, const ContigRegion& region);
};

// non-member methods

namespace detail {

template <typename InputIt>
void extend_tree(InputIt first, InputIt last, HaplotypeTree& tree, std::false_type)
{
    std::for_each(first, last, [&] (const auto& allele_or_haplotype) { tree.extend(allele_or_haplotype); });
}

template <typename InputIt>
void extend_tree(InputIt first, InputIt last, HaplotypeTree& tree, std::true_type)
{
    std::for_each(first, last, [&] (const auto& variant) {
        tree.extend(variant.ref_allele());
        tree.extend(variant.alt_allele());
    });
}

template <typename InputIt, typename A>
InputIt extend_tree_until(InputIt first, InputIt last, HaplotypeTree& tree,
                          const unsigned max_haplotypes, A, std::input_iterator_tag)
{
    if (first == last) return last;
    const auto it = std::find_if(first, last, [&] (const auto& allele) {
                                     tree.extend(allele);
                                     return tree.num_haplotypes() >= max_haplotypes;
                                 });
    if (tree.num_haplotypes() == max_haplotypes) {
        return std::next(it);
    } else {
        return it;
    }
}

inline unsigned max_log_haplotypes_after_extension(const std::size_t& tree_size, const unsigned num_new_alleles)
{
    return std::log2(tree_size) + num_new_alleles;
}

inline unsigned max_log_haplotypes_after_extension(const HaplotypeTree& tree, const unsigned num_new_alleles)
{
    return max_log_haplotypes_after_extension(tree.num_haplotypes(), num_new_alleles);
}

template <typename RandomIt, typename A>
RandomIt extend_tree_until(RandomIt first, RandomIt last, HaplotypeTree& tree, const std::size_t max_haplotypes,
                           A, std::random_access_iterator_tag)
{
    if (max_log_haplotypes_after_extension(tree, std::distance(first, last)) <= std::log2(max_haplotypes)) {
        extend_tree(first, last, tree, std::false_type {});
        return last;
    } else {
        return extend_tree_until(first, last, tree, max_haplotypes, A {}, std::input_iterator_tag {});
    }
}

template <typename InputIt, typename A>
InputIt extend_tree_until(InputIt first, InputIt last, HaplotypeTree& tree, const std::size_t max_haplotypes, A)
{
    return extend_tree_until(first, last, tree, max_haplotypes, A {},
                             typename std::iterator_traits<InputIt>::iterator_category {});
}

template <typename InputIt>
InputIt extend_tree_until(InputIt first, InputIt last, HaplotypeTree& tree, const std::size_t max_haplotypes,
                          Variant, std::input_iterator_tag)
{
    if (first == last) return last;
    const auto it = std::find_if(first, last, [&] (const auto& variant) {
                                     tree.extend(variant.ref_allele());
                                     tree.extend(variant.alt_allele());
                                     return tree.num_haplotypes() >= max_haplotypes;
                                 });
    if (tree.num_haplotypes() == max_haplotypes) {
        return std::next(it);
    } else {
        return it;
    }
}

template <typename RandomIt>
RandomIt extend_tree_until(RandomIt first, RandomIt last, HaplotypeTree& tree, const std::size_t max_haplotypes,
                           Variant, std::random_access_iterator_tag)
{
    if (max_log_haplotypes_after_extension(tree, 2 * std::distance(first, last)) <= std::log2(max_haplotypes)) {
        extend_tree(first, last, tree, Variant {});
        return last;
    } else {
        return extend_tree_until(first, last, tree, max_haplotypes, Variant {}, std::input_iterator_tag {});
    }
}

template <typename InputIt>
InputIt extend_tree_until(InputIt first, InputIt last, HaplotypeTree& tree, const std::size_t max_haplotypes, Variant)
{
    return extend_tree_until(first, last, tree, max_haplotypes, Variant {},
                             typename std::iterator_traits<InputIt>::iterator_category {});
}

template <typename T>
constexpr bool is_variant_or_allele_or_haplotype = std::is_same<T, ContigAllele>::value
                                                || std::is_same<T, Allele>::value
                                                || std::is_same<T, Variant>::value
                                                || std::is_same<T, Haplotype>::value;

} // namespace detail

template <typename InputIt>
void extend_tree(InputIt first, InputIt last, HaplotypeTree& tree)
{
    using MappableType = std::decay_t<typename std::iterator_traits<InputIt>::value_type>;
    static_assert(detail::is_variant_or_allele_or_haplotype<MappableType>, "not Allele or Variant or Haplotype");
    detail::extend_tree(first, last, tree, std::is_same<MappableType, Variant> {});
}

template <typename Container>
void extend_tree(const Container& elements, HaplotypeTree& tree)
{
    extend_tree(std::cbegin(elements), std::cend(elements), tree);
}

template <typename InputIt>
InputIt extend_tree_until(InputIt first, InputIt last, HaplotypeTree& tree, const std::size_t max_haplotypes)
{
    using MappableType = std::decay_t<typename std::iterator_traits<InputIt>::value_type>;
    static_assert(detail::is_variant_or_allele_or_haplotype<MappableType>, "not Allele or Variant or Haplotype");
    return detail::extend_tree_until(first, last, tree, max_haplotypes, MappableType {});
}

template <typename Container>
auto extend_tree_until(const Container& elements, HaplotypeTree& tree, const std::size_t max_haplotypes)
{
    return extend_tree_until(std::cbegin(elements), std::cend(elements), tree, max_haplotypes);
}

template <typename Container>
void prune_all(const Container& haplotypes, HaplotypeTree& tree)
{
    for (const auto& haplotype : haplotypes) {
        tree.prune_all(haplotype);
    }
}

template <typename Container>
void prune_unique(const Container& haplotypes, HaplotypeTree& tree)
{
    for (const auto& haplotype : haplotypes) {
        tree.prune_unique(haplotype);
    }
}

template <typename Container>
void splice(const Container& alleles, HaplotypeTree& tree)
{
    for (const auto& allele : alleles) {
        tree.splice(allele);
    }
}

template <typename InputIt>
auto generate_all_haplotypes(InputIt first, InputIt last, const ReferenceGenome& reference)
{
    using MappableType = std::decay_t<typename std::iterator_traits<InputIt>::value_type>;
    static_assert(detail::is_variant_or_allele_or_haplotype<MappableType>, "not Allele or Variant");
    if (first == last) return std::vector<Haplotype> {};
    HaplotypeTree tree {contig_name(*first), reference};
    extend_tree(first, last, tree);
    return tree.extract_haplotypes();
}

template <typename Container>
auto generate_all_haplotypes(const Container& elements, const ReferenceGenome& reference)
{
    return generate_all_haplotypes(std::cbegin(elements), std::cend(elements), reference);
}

HaplotypeTree::HaplotypeLength min_haplotype_length(const HaplotypeTree& tree);
HaplotypeTree::HaplotypeLength min_haplotype_length(const HaplotypeTree& tree, const GenomicRegion& region);
HaplotypeTree::HaplotypeLength max_haplotype_length(const HaplotypeTree& tree);
HaplotypeTree::HaplotypeLength max_haplotype_length(const HaplotypeTree& tree, const GenomicRegion& region);
std::pair<HaplotypeTree::HaplotypeLength, HaplotypeTree::HaplotypeLength>
minmax_haplotype_lengths(const HaplotypeTree& tree);
std::pair<HaplotypeTree::HaplotypeLength, HaplotypeTree::HaplotypeLength>
minmax_haplotype_lengths(const HaplotypeTree& tree, const GenomicRegion& region);

namespace debug {

void write_dot(const HaplotypeTree& tree, const boost::filesystem::path& dest);

} // namespace debug

} // namespace coretools
} // namespace octopus

#endif
