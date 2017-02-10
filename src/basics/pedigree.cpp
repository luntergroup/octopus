// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "pedigree.hpp"

#include <algorithm>
#include <iterator>

namespace octopus {

// public methods

void Pedigree::add_relationship(Member parent, Member child)
{
    const auto parent_vertex = get_member_vertex(parent);
    const auto child_vertex  = get_member_vertex(child);
    boost::add_edge(parent_vertex, child_vertex, tree_);
}

boost::optional<const SampleName&> Pedigree::mother_of(const SampleName& child) const
{
    return parent_of(child, Member::Sex::female);
}

boost::optional<const SampleName&> Pedigree::father_of(const SampleName& child) const
{
    return parent_of(child, Member::Sex::male);
}

std::vector<Pedigree::Member> Pedigree::children_of(const Member& parent) const
{
    const auto children = boost::adjacent_vertices(members_.at(parent.name), tree_);
    std::vector<Member> result {};
    result.reserve(std::distance(children.first, children.second));
    std::transform(children.first, children.second, std::back_inserter(result),
                   [this] (const Vertex& v) { return tree_[v]; });
    return result;
}

bool Pedigree::is_member(const SampleName& member) const noexcept
{
    return members_.count(member) == 1;
}

bool Pedigree::is_empty() const noexcept
{
    return members_.empty();
}

std::size_t Pedigree::size() const noexcept
{
    return members_.size();
}

void Pedigree::clear() noexcept
{
    tree_.clear();
    members_.clear();
}

// private methods

Pedigree::Vertex Pedigree::add_new_member(const Member& member)
{
    const auto result = boost::add_vertex(member, tree_);
    members_.insert(std::make_pair(member.name, result));
    return result;
}

Pedigree::Vertex Pedigree::get_member_vertex(const Member& member)
{
    if (!is_member(member.name)) {
        return add_new_member(member);
    } else {
        return members_.at(member.name);
    }
}

boost::optional<const SampleName&> Pedigree::parent_of(const SampleName& child, Member::Sex parent_gender) const
{
    const auto parents = boost::inv_adjacent_vertices(members_.at(child), tree_);
    const auto itr = std::find_if(parents.first, parents.second,
                                  [&] (const Vertex& v) { return tree_[v].gender == parent_gender; });
    if (itr == parents.second) {
        return boost::none;
    } else {
        return tree_[*itr].name;
    }
}

// non-member methods

boost::optional<Trio> make_trio(const SampleName& child, const Pedigree& pedigree)
{
    const auto mother = pedigree.mother_of(child);
    const auto father = pedigree.father_of(child);
    if (mother && father) {
        return Trio {{*mother}, {*father}, {child}};
    } else {
        return boost::none;
    }
}
    
} // namespace octopus