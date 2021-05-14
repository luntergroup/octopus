// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef pedigree_hpp
#define pedigree_hpp

#include <vector>
#include <unordered_map>
#include <cstddef>

#include <boost/graph/adjacency_list.hpp>
#include <boost/optional.hpp>

#include "config/common.hpp"
#include "trio.hpp"

namespace octopus {

class Pedigree
{
public:
    struct Member
    {
        enum class Sex { male, female, hermaphroditic };
        SampleName name;
        Sex sex = Sex::hermaphroditic;
    };
    
    Pedigree() = default;
    Pedigree(std::size_t pedigree_size_hint);
    
    Pedigree(const Pedigree&);
    Pedigree& operator=(const Pedigree&);
    Pedigree(Pedigree&&);
    Pedigree& operator=(Pedigree&&);
    
    ~Pedigree() = default;
    
    void add_founder(Member parent);
    void add_descendant(Member offspring, const SampleName& parent);
    void add_descendant(Member offspring, const SampleName& mother, const SampleName& father);
    
    bool is_member(const SampleName& member) const noexcept;
    bool is_founder(const SampleName& member) const;
    std::size_t num_parents(const SampleName& member) const;
    std::size_t num_offspring(const SampleName& member) const;

    boost::optional<const SampleName&> mother_of(const SampleName& child) const;
    boost::optional<const SampleName&> father_of(const SampleName& child) const;
    std::vector<Member> offspring_of(const Member& parent) const;
    
    bool is_empty() const noexcept;
    std::size_t size() const noexcept;
    void clear() noexcept;
    
private:
    using Tree = boost::adjacency_list<boost::listS, boost::listS, boost::bidirectionalS, Member, boost::no_property>;
    using Vertex = typename boost::graph_traits<Tree>::vertex_descriptor;
    using Edge   = typename boost::graph_traits<Tree>::edge_descriptor;
    
    Tree tree_;
    std::unordered_map<SampleName, Vertex> members_;
    
    Vertex add_new_member(Member member);
    Vertex get_vertex(const SampleName& member) const;
    boost::optional<const SampleName&> parent_of(const SampleName& child, Member::Sex parent_gender) const;
};

bool is_parent_of(const SampleName& parent, const SampleName& offspring, const Pedigree& pedigree);

bool is_trio(const std::vector<SampleName>& samples, const Pedigree& pedigree);

boost::optional<Trio> make_trio(const SampleName& child, const Pedigree& pedigree);

} // namespace octopus

#endif
