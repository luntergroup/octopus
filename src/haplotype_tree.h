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
    
    unsigned num_haplotypes() const;
    void extend_haplotypes(const Variant& a_variant);
    Haplotypes get_haplotypes(unsigned num_alleles);
    void prune_haplotypes(const Haplotypes& haplotypes);
    
private:
    struct AlleleNode
    {
        Allele the_allele;
    };
    
    using Tree = boost::adjacency_list<
        boost::listS, boost::listS, boost::undirectedS, AlleleNode, boost::no_property
    >;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree the_tree_;
    Vertex the_root_;
    std::list<Vertex> haplotype_branch_ends_;
    ReferenceGenome& the_reference_;
    std::unordered_map<Haplotype, Vertex> haplotype_to_branch_end_map_;
    unsigned haplotype_allele_length_;
    
    using BranchIterator = decltype(haplotype_branch_ends_)::const_iterator;
    
    Vertex get_previous_allele(Vertex allele) const;
    bool node_has_allele_branch(Vertex allele, const Allele& an_allele) const;
    BranchIterator extend_haplotype(BranchIterator haplotype_branch_end, const Variant& a_variant);
    BranchIterator extend_haplotype(BranchIterator haplotype_branch_end, const Allele& an_allele);
    bool is_haplotype_branch_end(Vertex haplotype_branch_end, const Haplotype& haplotype) const;
    Vertex find_haplotype_branch_end(const Haplotype& haplotype) const;
    void prune_haplotype(const Haplotype& haplotype);
    void prune_haplotype_from_branch_end(Vertex haplotype_branch_end);
};

#endif /* defined(__Octopus__haplotype_tree__) */
