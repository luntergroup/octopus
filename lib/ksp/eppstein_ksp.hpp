// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef EPPSTEIN_KSP_HPP
#define EPPSTEIN_KSP_HPP

#include <boost/graph/graph_traits.hpp>

template <typename VertexListGraph, typename WeightMap, typename IndexMap, typename ResultMap>
void
eppstein_ksp(const VertexListGraph& g,
             typename boost::graph_traits<VertexListGraph>::vertex_descriptor s,
             typename boost::graph_traits<VertexListGraph>::vertex_descriptor t,
             WeightMap weight, IndexMap index_map, unsigned k, ResultMap& result)
{
    
}

template <typename VertexListGraph, typename WeightMap, typename IndexMap>
auto
eppstein_ksp(const VertexListGraph& g,
             typename boost::graph_traits<VertexListGraph>::vertex_descriptor s,
             typename boost::graph_traits<VertexListGraph>::vertex_descriptor t,
             WeightMap weight, IndexMap index, unsigned k)
{
    using WeightType = typename WeightMap::mapped_type;
    using EdgeType = typename boost::graph_traits<VertexListGraph>::edge_descriptor;
    std::map<WeightType, std::deque<EdgeType>> result {};
    eppstein_ksp(g, s, t, std::move(weight), std::move(index), k, result);
    return result;
}

#endif //EPPSTEIN_KSP_HPP
