//
//  haplotype_tree.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_tree.h"

#include <stdexcept>
#include <iterator>  // std::next, std::tie, std::cbegin etc
#include <algorithm> // std::find, std::find_if

#include "reference_genome.h"
#include "variant.h"
#include "region_utils.h"

#include <iostream> // TEST

HaplotypeTree::HaplotypeTree(ReferenceGenome& the_reference)
:
the_reference_ {the_reference},
the_tree_ {},
the_root_ {boost::add_vertex(the_tree_)},
haplotype_leafs_ {the_root_},
haplotype_leaf_cache_ {},
haplotype_allele_length_ {}
{}

unsigned HaplotypeTree::num_haplotypes() const
{
    return static_cast<unsigned>(haplotype_leafs_.size());
}

void HaplotypeTree::extend_haplotypes(const Allele& an_allele)
{
    auto branch_end_it = std::cbegin(haplotype_leafs_);
    
    while (branch_end_it != std::cend(haplotype_leafs_)) {
        branch_end_it = extend_haplotype(branch_end_it, an_allele);
        ++branch_end_it;
    }
}

HaplotypeTree::Haplotypes HaplotypeTree::get_haplotypes(const GenomicRegion& a_region)
{
    haplotype_leaf_cache_.clear(); // to prevent the cache becoming out of date
    
    Haplotypes result {};
    
    for (auto leaf : haplotype_leafs_) {
        auto haplotype = get_haplotype(leaf, a_region);
        
        // recently retreived haplotypes are added to the cache as it is likely these
        // are the haplotypes that will be pruned next
        haplotype_leaf_cache_.emplace(haplotype, leaf);
        
        result.emplace_back(std::move(haplotype));
    }
    
    return result;
}

void HaplotypeTree::prune_haplotype(const Haplotype& haplotype)
{
    auto max_occurences = num_haplotypes();
    
    Vertex new_haplotype_end;
    bool new_end_is_leaf;
    
    auto leaf_it = std::cbegin(haplotype_leafs_);
    
    for (unsigned i {}; i < max_occurences; ++i) {
        
        if (haplotype_leaf_cache_.count(haplotype) > 0) {
            leaf_it = std::find(leaf_it, std::cend(haplotype_leafs_), haplotype_leaf_cache_.at(haplotype));
            haplotype_leaf_cache_.erase(haplotype);
        } else {
            leaf_it = find_haplotype_leaf(leaf_it, std::cend(haplotype_leafs_), haplotype);
        }
        
        if (leaf_it == std::cend(haplotype_leafs_)) {
            return;
        }
        
        std::tie(new_haplotype_end, new_end_is_leaf) = prune_branch(*leaf_it, haplotype.get_region());
        
        leaf_it = haplotype_leafs_.erase(leaf_it);
        
        if (new_end_is_leaf) {
            leaf_it = haplotype_leafs_.insert(leaf_it, new_haplotype_end);
        }
    }
}

HaplotypeTree::Vertex HaplotypeTree::get_previous_allele(Vertex allele) const
{
    return *boost::inv_adjacent_vertices(allele, the_tree_).first;
}

bool HaplotypeTree::allele_exists(Vertex allele, const Allele& an_allele) const
{
    auto vertex_range = boost::adjacent_vertices(allele, the_tree_);
    return std::any_of(vertex_range.first, vertex_range.second, [this, &an_allele] (const auto& vertex) {
        return the_tree_[vertex].the_allele == an_allele;
    });
}

HaplotypeTree::LeafIterator HaplotypeTree::extend_haplotype(LeafIterator haplotype_leaf, const Variant& a_variant)
{
    haplotype_leaf = extend_haplotype(haplotype_leaf, a_variant.get_reference_allele());
    return extend_haplotype(haplotype_leaf, a_variant.get_alternative_allele());
}

HaplotypeTree::LeafIterator HaplotypeTree::extend_haplotype(LeafIterator haplotype_leaf,
                                                              const Allele& the_new_allele)
{
    const auto& branch_allele = the_tree_[*haplotype_leaf].the_allele;
    
    if (*haplotype_leaf != the_root_ && is_before(the_new_allele, branch_allele)) {
        return haplotype_leaf;
    } else if (*haplotype_leaf == the_root_ || !overlaps(branch_allele, the_new_allele)) {
        auto new_branch = boost::add_vertex(the_tree_);
        the_tree_[new_branch].the_allele = the_new_allele;
        boost::add_edge(*haplotype_leaf, new_branch, the_tree_);
        haplotype_leaf = haplotype_leafs_.erase(haplotype_leaf);
        haplotype_leaf = haplotype_leafs_.insert(haplotype_leaf, new_branch);
    } else {
        if (begins_equal(branch_allele, the_new_allele)) {
            Vertex previous_allele = get_previous_allele(*haplotype_leaf);
            if (!allele_exists(previous_allele, the_new_allele)) {
                auto new_branch = boost::add_vertex(the_tree_);
                the_tree_[new_branch].the_allele = the_new_allele;
                boost::add_edge(previous_allele, new_branch, the_tree_);
                haplotype_leafs_.insert(haplotype_leaf, new_branch);
            }
        }
    }
    
    return haplotype_leaf;
}

Haplotype HaplotypeTree::get_haplotype(Vertex haplotype_end, const GenomicRegion& a_region)
{
    Haplotype result {the_reference_};
    
    while (haplotype_end != the_root_ && is_after(the_tree_[haplotype_end].the_allele, a_region)) {
        haplotype_end = get_previous_allele(haplotype_end);
    }
    
    while (haplotype_end != the_root_ && overlaps(the_tree_[haplotype_end].the_allele, a_region)) {
        result.push_front(the_tree_[haplotype_end].the_allele);
        haplotype_end = get_previous_allele(haplotype_end);
    }
    
    return result;
}

bool HaplotypeTree::is_branch_the_haplotype(Vertex haplotype_end, const Haplotype& haplotype) const
{
    if (!overlaps(haplotype, the_tree_[haplotype_end].the_allele)) {
        return false;
    }
    
    while (haplotype_end != the_root_ && !is_before(the_tree_[haplotype_end].the_allele, haplotype)) {
        if (!haplotype.contains(the_tree_[haplotype_end].the_allele)) return false;
        haplotype_end = get_previous_allele(haplotype_end);
    }
    
    return true;
}

HaplotypeTree::LeafIterator HaplotypeTree::find_haplotype_leaf(LeafIterator first, LeafIterator last,
                                                               const Haplotype& haplotype) const
{
    return std::find_if(first, last, [this, &haplotype] (Vertex leaf) {
        return is_branch_the_haplotype(leaf, haplotype);
    });
}

std::pair<HaplotypeTree::Vertex, bool> HaplotypeTree::prune_branch(Vertex leaf, const GenomicRegion& a_region)
{
    while (leaf != the_root_) {
        if (boost::out_degree(leaf, the_tree_) > 0) {
            return {leaf, false};
        } else if (begins_before(the_tree_[leaf].the_allele, a_region)) {
            return {leaf, true};
        }
        
        auto new_leaf = get_previous_allele(leaf);
        boost::remove_edge(new_leaf, leaf, the_tree_);
        boost::remove_vertex(leaf, the_tree_);
        leaf = new_leaf;
    }
    
    return {leaf, true};
}
