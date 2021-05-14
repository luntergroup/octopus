// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "pedigree.hpp"

#include <utility>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <stdexcept>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/copy.hpp>

namespace octopus {

// public methods

Pedigree::Pedigree(std::size_t pedigree_size_hint)
{
    members_.reserve(pedigree_size_hint);
}

namespace {

template <typename Graph>
auto copy_graph(const Graph& src, Graph& dst)
{
    assert(boost::num_vertices(dst) == 0);
    using Vertex = typename boost::graph_traits<Graph>::vertex_descriptor;
    std::unordered_map<Vertex, std::size_t> index_map {};
    index_map.reserve(boost::num_vertices(src));
    const auto p = boost::vertices(src);
    std::size_t i {0};
    std::for_each(p.first, p.second, [&i, &index_map] (Vertex v) { index_map.emplace(v, i++); });
    std::unordered_map<Vertex, Vertex> vertex_copy_map {};
    vertex_copy_map.reserve(boost::num_vertices(src));
    boost::copy_graph(src, dst,
                      boost::vertex_index_map(boost::make_assoc_property_map(index_map))
                      .orig_to_copy(boost::make_assoc_property_map(vertex_copy_map)));
    assert(vertex_copy_map.size() == boost::num_vertices(src));
    return vertex_copy_map;
}

template <typename Map1, typename Map2>
void copy_members(const Map1& src, Map1& dst, const Map2& vertex_copy_map)
{
    assert(dst.empty());
    dst.reserve(src.size());
    std::transform(std::cbegin(src), std::cend(src), std::inserter(dst, std::begin(dst)),
                   [&vertex_copy_map] (const auto& p) {
                       return std::make_pair(p.first, vertex_copy_map.at(p.second));
                   });
}

} // namespace

Pedigree::Pedigree(const Pedigree& other)
: tree_ {}
{
    const auto vertex_copy_map = copy_graph(other.tree_, tree_);
    copy_members(other.members_, members_, vertex_copy_map);
}

Pedigree& Pedigree::operator=(const Pedigree& other)
{
    if (&other == this) return *this;
    tree_.clear();
    members_.clear();
    const auto vertex_copy_map = copy_graph(other.tree_, tree_);
    copy_members(other.members_, members_, vertex_copy_map);
    return *this;
}

Pedigree::Pedigree(Pedigree&& other) : Pedigree {other} {}

Pedigree& Pedigree::operator=(Pedigree&& other)
{
    *this = other;
    return *this;
}

void Pedigree::add_founder(Member parent)
{
     add_new_member(std::move(parent));
}

void Pedigree::add_descendant(Member offspring, const SampleName& parent)
{
    const auto offspring_vertex = add_new_member(std::move(offspring));
    try {
        const auto parent_vertex = get_vertex(parent);
        boost::add_edge(parent_vertex, offspring_vertex, tree_);
    } catch (const std::out_of_range&) {
        throw std::runtime_error {"Parent not in pedigree"};
    }
}

void Pedigree::add_descendant(Member offspring, const SampleName& mother, const SampleName& father)
{
    const auto offspring_vertex = add_new_member(std::move(offspring));
    try {
        const auto mother_vertex = get_vertex(mother);
        const auto father_vertex = get_vertex(father);
        boost::add_edge(mother_vertex, offspring_vertex, tree_);
        boost::add_edge(father_vertex, offspring_vertex, tree_);
    } catch (const std::out_of_range&) {
        throw std::runtime_error {"Parent not in pedigree"};
    }
}

bool Pedigree::is_member(const SampleName& member) const noexcept
{
    return members_.count(member) == 1;
}

bool Pedigree::is_founder(const SampleName& member) const
{
    return num_parents(member) == 0;
}

std::size_t Pedigree::num_parents(const SampleName& member) const
{
    const auto parents = boost::inv_adjacent_vertices(get_vertex(member), tree_);
    return std::distance(parents.first, parents.second);
}

std::size_t Pedigree::num_offspring(const SampleName& member) const
{
    const auto offspring = boost::adjacent_vertices(get_vertex(member), tree_);
    return std::distance(offspring.first, offspring.second);
}

boost::optional<const SampleName&> Pedigree::mother_of(const SampleName& child) const
{
    return parent_of(child, Member::Sex::female);
}

boost::optional<const SampleName&> Pedigree::father_of(const SampleName& child) const
{
    return parent_of(child, Member::Sex::male);
}

std::vector<Pedigree::Member> Pedigree::offspring_of(const Member& parent) const
{
    const auto offspring = boost::adjacent_vertices(get_vertex(parent.name), tree_);
    std::vector<Member> result {};
    result.reserve(std::distance(offspring.first, offspring.second));
    std::transform(offspring.first, offspring.second, std::back_inserter(result),
                   [this] (const Vertex& v) { return tree_[v]; });
    return result;
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

Pedigree::Vertex Pedigree::add_new_member(Member member)
{
    const auto result = boost::add_vertex(member, tree_);
    members_.insert(std::make_pair(std::move(member.name), result));
    return result;
}

Pedigree::Vertex Pedigree::get_vertex(const SampleName& member) const
{
    return members_.at(member);
}

boost::optional<const SampleName&> Pedigree::parent_of(const SampleName& child, const Member::Sex parent_sex) const
{
    const auto parents = boost::inv_adjacent_vertices(get_vertex(child), tree_);
    const auto itr = std::find_if(parents.first, parents.second, [&] (Vertex v) { return tree_[v].sex == parent_sex; });
    if (itr == parents.second) {
        return boost::none;
    } else {
        return tree_[*itr].name;
    }
}

// non-member methods

bool is_parent_of(const SampleName& parent, const SampleName& offspring, const Pedigree& pedigree)
{
    const auto mother = pedigree.mother_of(offspring);
    if (mother && *mother == parent) return true;
    const auto father = pedigree.father_of(offspring);
    return father && *father == parent;
}

bool all_members(const std::vector<SampleName>& samples, const Pedigree& pedigree)
{
    return std::all_of(std::cbegin(samples), std::cend(samples),
                       [&] (const auto& sample) { return pedigree.is_member(sample); });
}

bool is_trio(const std::vector<SampleName>& samples, const Pedigree& pedigree)
{
    if (samples.size() == 3 && all_members(samples, pedigree)) {
        if (is_parent_of(samples[0], samples[2], pedigree)) {
            return is_parent_of(samples[1], samples[2], pedigree);
        } else if (is_parent_of(samples[0], samples[1], pedigree)) {
            return is_parent_of(samples[2], samples[1], pedigree);
        } else {
            return is_parent_of(samples[1], samples[0], pedigree)
                   && is_parent_of(samples[2], samples[0], pedigree);
        }
    } else {
        return false;
    }
}

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