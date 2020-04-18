//
// Created by fkuner on 2020/4/7.
//

#ifndef TUX_INGRESS_EDGE_DECISION_HPP
#define TUX_INGRESS_EDGE_DECISION_HPP

#include <tux/graph/distributed_graph.hpp>
#include <tux/graph/graph_basic_types.hpp>

#include <graphlab/graph/graph_hash.hpp>
#include <graphlab/rpc/distributed_event_log.hpp>
#include <graphlab/util/dense_bitset.hpp>
#include <boost/random/uniform_int_distribution.hpp>

namespace tux{
    template <typename FVertexType, typename SVertexType, typename EdgeType>
    class distributed_graph;

    template <typename FVertexType, typename SVertexType, typename EdgeType>
    class ingress_edge_decision{

    };
}

#endif //TUX_INGRESS_EDGE_DECISION_HPP
