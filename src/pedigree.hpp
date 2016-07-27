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
#include <cstddef>

#include <boost/graph/adjacency_list.hpp>

#include "common.hpp"

class Pedigree
{
public:
    using Member = Octopus::SampleIdType;
    
    Pedigree() = default;
    
    Pedigree(const Pedigree&)            = default;
    Pedigree& operator=(const Pedigree&) = default;
    Pedigree(Pedigree&&)                 = default;
    Pedigree& operator=(Pedigree&&)      = default;
    
    ~Pedigree() = default;
    
    void add_relationship(Member parent, Member child);
    void clear();
    std::size_t size() const;
    
    bool is_child(Member child, Member parent) const;
    bool is_descendant(Member descendant, Member ancestor) const;
    unsigned num_children(Member parent) const;
    std::vector<Member> get_children(Member parent) const;
    std::vector<Member> get_descendants(Member parent) const;
    
private:
    using Tree = boost::adjacency_list<
                    boost::listS, boost::listS, boost::bidirectionalS, Member, boost::no_property
                >;
    
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree tree_;
    
    std::vector<Member> roots_;
    std::vector<Member> leafs_;
};

#endif /* pedigree_hpp */
