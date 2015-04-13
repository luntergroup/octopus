//
//  haplotype_tree.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_tree.h"

#include <stdexcept>
#include <iterator> // std::next, std::cbegin etc

#include "reference_genome.h"
#include "variant.h"
#include "region_utils.h"

#include <iostream> // TEST

HaplotypeTree::HaplotypeTree(ReferenceGenome& the_reference)
:
the_reference_ {the_reference},
the_tree_ {},
the_root_ {boost::add_vertex(the_tree_)},
haplotype_branch_ends_ {the_root_},
haplotype_to_branch_end_map_ {},
haplotype_allele_length_ {}
{}

unsigned HaplotypeTree::num_haplotypes() const
{
    return static_cast<unsigned>(haplotype_branch_ends_.size());
}

void HaplotypeTree::extend_haplotypes(const Variant& a_variant)
{
    auto branch_end_it      = std::cbegin(haplotype_branch_ends_);
    auto last_branch_end_it = std::cend(haplotype_branch_ends_);
    
    while (branch_end_it != last_branch_end_it) {
        branch_end_it = extend_haplotype(branch_end_it, a_variant);
    }
}

HaplotypeTree::Haplotypes HaplotypeTree::get_haplotypes(unsigned num_alleles)
{
    Haplotypes result {};
    
    for (auto haplotype_end : haplotype_branch_ends_) {
        Haplotype a_haplotype {the_reference_};
        
        for (size_t i {0}; i < num_alleles; ++i) {
            a_haplotype.push_front(the_tree_[haplotype_end].the_allele);
            haplotype_end = get_previous_allele(haplotype_end);
        }
        
        haplotype_to_branch_end_map_[a_haplotype] = haplotype_end;
        result.emplace_back(std::move(a_haplotype));
    }
    
    return result;
}

void HaplotypeTree::prune_haplotypes(const Haplotypes& haplotypes)
{
    for (const auto& haplotype : haplotypes) {
        prune_haplotype(haplotype);
    }
}

HaplotypeTree::Vertex HaplotypeTree::get_previous_allele(Vertex allele) const
{
    return *boost::inv_adjacent_vertices(allele, the_tree_).first;
}

bool HaplotypeTree::node_has_allele_branch(Vertex allele, const Allele& an_allele) const
{
    auto vertex_range = boost::adjacent_vertices(allele, the_tree_);
    return std::any_of(vertex_range.first, vertex_range.second, [this, &an_allele] (const auto& vertex) {
        return the_tree_[vertex].the_allele == an_allele;
    });
}

HaplotypeTree::BranchIterator HaplotypeTree::extend_haplotype(BranchIterator haplotype_branch_end, const Variant& a_variant)
{
    haplotype_branch_end = extend_haplotype(haplotype_branch_end, a_variant.get_reference_allele());
    return extend_haplotype(haplotype_branch_end, a_variant.get_alternative_allele());
}

HaplotypeTree::BranchIterator HaplotypeTree::extend_haplotype(BranchIterator haplotype_branch_end, const Allele& an_allele)
{
    if (!overlaps(the_tree_[*haplotype_branch_end].the_allele, an_allele)) {
        auto new_branch = boost::add_vertex(the_tree_);
        the_tree_[new_branch].the_allele = an_allele;
        boost::add_edge(*haplotype_branch_end, new_branch, the_tree_);
        haplotype_branch_end = haplotype_branch_ends_.erase(haplotype_branch_end);
        return haplotype_branch_ends_.insert(haplotype_branch_end, new_branch);
    } else {
        if (contains(an_allele, the_tree_[*haplotype_branch_end].the_allele)) {
            std::cout << "CONTAINS: " << an_allele << std::endl;
        } else {
            std::cout << "OVERLAP: " << an_allele << std::endl;
            
            Vertex previous_allele = get_previous_allele(*haplotype_branch_end);
//            while (previous_allele != the_root_ && ) {
//                
//            }
            
            if (begins_equal(the_tree_[previous_allele].the_allele, an_allele)) {
                if (!node_has_allele_branch(previous_allele, an_allele)) {
                    auto new_branch = boost::add_vertex(the_tree_);
                    the_tree_[new_branch].the_allele = an_allele;
                    boost::add_edge(*haplotype_branch_end, new_branch, the_tree_);
                    return haplotype_branch_ends_.insert(haplotype_branch_end, new_branch);
                }
            } else {
                
            }
        }
        
        return haplotype_branch_end;
    }
}

bool HaplotypeTree::is_haplotype_branch_end(Vertex haplotype_branch_end, const Haplotype& haplotype) const
{
    auto haplotype_region = haplotype.get_region();
    auto allele_vertex = haplotype_branch_end;
    
    while (allele_vertex != the_root_ && !is_before(the_tree_[allele_vertex].the_allele, haplotype_region)) {
        if (!haplotype.contains(the_tree_[allele_vertex].the_allele)) return false;
    }
    
    return true;
}

HaplotypeTree::Vertex HaplotypeTree::find_haplotype_branch_end(const Haplotype& haplotype) const
{
    if (haplotype_to_branch_end_map_.count(haplotype) > 0) {
       return haplotype_to_branch_end_map_.at(haplotype);
    } else {
        for (const auto& haplotype_branch_end : haplotype_branch_ends_) {
            if (is_haplotype_branch_end(haplotype_branch_end, haplotype)) {
                return haplotype_branch_end;
            }
        }
    }
    throw std::runtime_error {"Haplotype not in tree"};
}

void HaplotypeTree::prune_haplotype(const Haplotype& haplotype)
{
    prune_haplotype_from_branch_end(find_haplotype_branch_end(haplotype));
}

void HaplotypeTree::prune_haplotype_from_branch_end(Vertex haplotype_branch_end)
{
    while (boost::out_degree(haplotype_branch_end, the_tree_) == 1) {
        auto new_haplotype_end = *boost::inv_adjacent_vertices(haplotype_branch_end, the_tree_).first;
        boost::remove_edge(new_haplotype_end, haplotype_branch_end, the_tree_);
        boost::remove_vertex(haplotype_branch_end, the_tree_);
        haplotype_branch_end = new_haplotype_end;
    }
}
