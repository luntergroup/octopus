// =======================================================================
// Copyright 2015 by Ireneusz Szcześniak
// Author: Ireneusz Szcześniak <www.irkos.org>
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
// =======================================================================

// =======================================================================
// This is the implementation of the Yen algorithm:
//
// Jin Y. Yen, Finding the k shortest loopless paths in a network,
// Management Science, vol. 17, no. 11, July 1971, pages 712-716
//
// Note 1: in the article there are Q^k Q_k Q^k_k used.  As far as I
// can tell they have the same meaning.  Q simply denotes the number
// of the node previous to the last node in a path (which always
// equals to the number of nodes in the path - 1), while (Q) denotes
// the node number Q in a path.  For instance, in the path a-b-c-d,
// the number of the node previous to the last is Q = 3, while (Q) =
// c.  Q^k refers to the k-th shortest path, and simply means the
// number of the node, which is previous to the last in the k-shortest
// path.
// =======================================================================

#ifndef BOOST_GRAPH_YEN_KSP
#define BOOST_GRAPH_YEN_KSP

#include <list>
#include <set>
#include <utility>

#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/optional.hpp>

#include "custom_dijkstra_call.hpp"

namespace boost {

  template <typename Graph, typename WeightMap, typename IndexMap>
  std::list<std::pair<typename WeightMap::value_type,
                      std::list<typename Graph::edge_descriptor>>>
  yen_ksp(const Graph& g,
          typename Graph::vertex_descriptor s,
          typename Graph::vertex_descriptor t,
          WeightMap wm, IndexMap im, optional<unsigned> K)
  {
    typedef typename Graph::vertex_descriptor vertex_descriptor;
    typedef typename Graph::edge_descriptor edge_descriptor;
    typedef typename WeightMap::value_type weight_type;
    typedef std::set<edge_descriptor> es_type;
    typedef std::set<vertex_descriptor> vs_type;
    typedef filtered_graph<Graph, is_not_in_subset<es_type>,
                           is_not_in_subset<vs_type>> fg_type;
    typedef std::list<edge_descriptor> path_type;
    typedef std::pair<weight_type, path_type> kr_type;

    // The shortest paths - these we return.
    std::list<kr_type> A;

    // An empty result if the source and destination are the same.
    if (s == t)
      return A;

    // The tentative paths - these are candidate paths.  It's a set,
    // because we want to make sure that a given result can show up in
    // the set of tentative results only once.  The problem is that
    // the algorithm can find the same tentative path many times.
    std::set<kr_type> B;

    // Try to find the (optional) shortest path.
    optional<kr_type> osp = custom_dijkstra_call(g, s, t, wm, im);

    // We quit if there was no shortest path.
    if (!osp)
      return A;
    
    // The first shortest path found becomes our first solution.
    A.push_back(std::move(osp.get()));

    // In each iteration we produce the next k-th shortest path.
    for (int k = 2; !K || k <= static_cast<int>(K.get()); ++k)
      {
        // The previous shortest result and path.
        const auto &psr = A.back();
        const auto &psp = psr.second;

        // The set of excluded edges.  It's a set, because the
        // algorithm can try to exclude an edge many times.
        es_type exe;
        // The set of excluded vertexes.
        vs_type exv;

        // The edge predicate.
        is_not_in_subset<es_type> ep(exe);
        // The vertex predicate.
        is_not_in_subset<vs_type> vp(exv);

        // The filtered graph.
        fg_type fg(g, ep, vp);

        // The root result: the cost and the root path.
        kr_type rr;
        // The root path.
        const path_type &rp = rr.second;

        // Use the previous shortest path to get tentative paths.  We
        // can go ahead with the loop without checking any condition
        // (the condition in the for-statement is true): the path
        // found must have at least one link, because s != t.
        for(typename path_type::const_iterator i = psp.begin(); true;)
          {
            // An edge of the previous shortest path.
            const edge_descriptor &edge = *i;

            // The spur vertex - we try to deviate at this node.
            const vertex_descriptor &sv = source(edge, g);

            // Iterate over all previous shortest paths.
            // An iteration examines j-th shortest path.
            for(const auto &jr: A)
              {
                // The j-th shortest path.
                const path_type &jp = jr.second;

                // Let's prepare for the comparison.
                typename path_type::const_iterator jpi = jp.begin();
                typename path_type::const_iterator rpi = rp.begin();

                // Iterate as long as the edges are equal.
                while(jpi != jp.end() && rpi != rp.end() && *jpi == *rpi)
                  ++jpi, ++rpi;

                // Make sure we didn't reach the end of jp.  If we
                // did, there is no next edge in jp, which we could
                // exclude.  Also, make sure we reached the end of rp,
                // i.e., the jp begins with the complete rp, and not a
                // head of rp.
                if (jpi != jp.end() && rpi == rp.end())
                  exe.insert(*jpi);
              }

            // Optional spur result.
            optional<kr_type> osr = custom_dijkstra_call(fg, sv, t, wm, im);

            if (osr)
              {
                // The tentative result.
                kr_type tr = std::move(osr.get());
                tr.first += rr.first;
                tr.second.insert(tr.second.begin(), rp.begin(), rp.end());
                B.insert(std::move(tr));
              }

            // We have the condition to break the look here, and not
            // at the beginning, because we don't want to execute the
            // remainer of the loop in vain.
            if (++i == psp.end())
              break;

            // Remove the vertex that in this iteration is the spur,
            // but in the next iteration it's going to be a vertex
            // that should not be considered in the search for a spur
            // path.
            exv.insert(sv);

            // Add the edge to the back of the root result.
            rr.first += get(wm, edge);
            rr.second.push_back(edge);
          }

        // Stop searching when there are no tentative paths.
        if (B.empty())
          break;

        // Take the shortest tentative path and make it the next
        // shortest path.
        A.push_back(std::move(*B.begin()));
        B.erase(B.begin());
      }
    
    return A;
  }

  template <typename Graph>
  std::list<std::pair<typename property_map<Graph, edge_weight_t>::value_type,
                      std::list<typename Graph::edge_descriptor>>>
  yen_ksp(Graph& g,
          typename Graph::vertex_descriptor s,
          typename Graph::vertex_descriptor t,
          optional<unsigned> K = optional<unsigned>())
  {
    return yen_ksp(g, s, t, get(edge_weight_t(), g),
                   get(vertex_index_t(), g), K);
  }

} // boost

#endif /* BOOST_GRAPH_YEN_KSP */
