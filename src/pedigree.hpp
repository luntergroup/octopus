//
//  pedigree.hpp
//  Octopus
//
//  Created by Daniel Cooke on 23/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef pedigree_hpp
#define pedigree_hpp

#include <vector>

#include <boost/graph/adjacency_list.hpp>

#include "common.hpp"

class Pedigree
{
public:
    using Member = Octopus::SampleIdType;
    
    Pedigree() = delete;
    Pedigree(Member root_mother, Member root_father);
    ~Pedigree() = default;
    
    Pedigree(const Pedigree&)            = default;
    Pedigree& operator=(const Pedigree&) = default;
    Pedigree(Pedigree&&)                 = default;
    Pedigree& operator=(Pedigree&&)      = default;
    
    void add_relationship(Member mother, Member father, Member child);
    
    bool is_child(Member child, Member parent) const;
    bool is_descendant(Member descendant, Member ancestor) const;
    unsigned num_children(Member parent) const;
    std::vector<Member> get_children(Member parent) const;
    std::vector<Member> get_descendants(Member parent) const;
    
private:
    using Tree   = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, Member, boost::no_property>;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree tree_;
    
    
};

#endif /* pedigree_hpp */
