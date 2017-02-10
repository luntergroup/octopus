// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef pedigree_hpp
#define pedigree_hpp

#include <vector>
#include <unordered_map>
#include <cstddef>

#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>

#include "config/common.hpp"
#import "trio.hpp"

namespace octopus {

class Pedigree
{
public:
    struct Member
    {
        enum class Sex { male, female };
        SampleName name;
        Sex gender;
    };
    
    Pedigree() = default;
    
    Pedigree(const Pedigree&)            = default;
    Pedigree& operator=(const Pedigree&) = default;
    Pedigree(Pedigree&&)                 = default;
    Pedigree& operator=(Pedigree&&)      = default;
    
    ~Pedigree() = default;
    
    void add_relationship(Member parent, Member child);
    
    boost::optional<const SampleName&> mother_of(const SampleName& child) const;
    boost::optional<const SampleName&> father_of(const SampleName& child) const;
    std::vector<Member> children_of(const Member& parent) const;
    
    bool is_member(const SampleName& member) const noexcept;
    
    bool is_empty() const noexcept;
    std::size_t size() const noexcept;
    void clear() noexcept;
    
private:
    using Tree = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, Member, boost::no_property>;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree tree_;
    std::unordered_map<SampleName, Vertex> members_;
    
    Vertex add_new_member(const Member& member);
    Vertex get_member_vertex(const Member& member);
    boost::optional<const SampleName&> parent_of(const SampleName& child, Member::Sex parent_gender) const;
};

boost::optional<Trio> make_trio(const SampleName& child, const Pedigree& pedigree);

} // namespace octopus

#endif
