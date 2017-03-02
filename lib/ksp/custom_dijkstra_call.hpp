// =======================================================================
// Copyright 2015 by Ireneusz Szcześniak
// Authors: Ireneusz Szcześniak
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
// =======================================================================

// =======================================================================
// The custom Dijkstra call, which returns the optional path as a list
// of edges along with the cost of the path.  The search is stopped
// when the destination node is reached.
// =======================================================================

#ifndef BOOST_GRAPH_CUSTOM_DIJKSTRA_CALL
#define BOOST_GRAPH_CUSTOM_DIJKSTRA_CALL

#include <list>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/visitors.hpp>
#include <boost/optional.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/utility/value_init.hpp>

namespace boost {

  // ========================================================================
  // Finish the search when a given node is examined, i.e. when the
  // shortest path to that node is found.
  // ========================================================================

  // The type of the exception thrown by the cdc_visitor.
  struct cdc_exception {};

  template <class Graph>
  struct cdc_visitor
  {
    typedef typename Graph::vertex_descriptor vertex_descriptor;
    typedef on_examine_vertex event_filter;
    cdc_visitor(vertex_descriptor t): m_t(t) {}
    void operator()(vertex_descriptor v, const Graph& g) {
      if (v == m_t)
        throw cdc_exception();
    }
    vertex_descriptor m_t;
  };

  // =======================================================================
  // The function.
  // =======================================================================

  template <typename Graph, typename WeightMap, typename IndexMap>
  optional<std::pair<typename WeightMap::value_type,
                     std::list<typename Graph::edge_descriptor>>>
  custom_dijkstra_call(const Graph &g,
                       typename Graph::vertex_descriptor s,
                       typename Graph::vertex_descriptor t,
                       WeightMap wm, IndexMap im)
  {
    typedef typename Graph::vertex_descriptor vertex_descriptor;
    typedef typename Graph::edge_descriptor edge_descriptor;
    typedef typename std::list<typename Graph::edge_descriptor> path_type;
    typedef typename WeightMap::value_type weight_type;
    typedef typename std::pair<weight_type, path_type> kr_type;

    optional<kr_type> result;

    if (s == t)
      result = std::make_pair(0, path_type());
    else
      {
        std::vector<edge_descriptor> pred_vec(num_vertices(g));
        auto pred = make_iterator_property_map(pred_vec.begin(), im);
        auto rep = record_edge_predecessors(pred, on_edge_relaxed());
        auto qat = cdc_visitor<Graph>(t);
        auto dv = make_dijkstra_visitor(std::make_pair(rep, qat));

        try
          {
            dijkstra_shortest_paths(g, s,
                                    weight_map(wm).vertex_index_map(im).
                                    visitor(dv));
          }
        catch (cdc_exception) {}

        // Was the solution found?
        if (pred[t] != edge_descriptor())
          {
            // The cost of the shortest path.
            value_initialized<weight_type> cost;
            // The path found.
            path_type p;

            // Trace the solution to the source.
            vertex_descriptor c = t;
            while (c != s)
              {
                const edge_descriptor &e = pred[c];
                // Build the path.
                p.push_front(e);
                // Calculate the cost of the path.
                cost += get(wm, e);
                // Find the predecessing vertex.
                c = source(e, g);
              }

            result = std::make_pair(cost, p);
          }
      }

    return result;
  }

} // boost

#endif /* BOOST_GRAPH_CUSTOM_DIJKSTRA_CALL */
