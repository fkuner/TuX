//
// Created by fkuner on 2020/3/20.
//

#ifndef TUX_GRAPH_BASIC_TYPES_HPP
#define TUX_GRAPH_BASIC_TYPES_HPP

#include <stdint.h>


namespace tux{

#ifdef USE_VID32
    /// Identifier type of a vertex which is globally consistent. Guaranteed to be integral
    typedef uint32_t vertex_id_type;
#else
    typedef uint64_t vertex_id_type;
#endif

    /// Identifier type of a vertex which is only locally consistent. Guaranteed to be integral
    typedef vertex_id_type lvid_type;

    /**
     * Identifier type of an edge which is only locally
     * consistent. Guaranteed to be integral and consecutive.
     */
    typedef lvid_type edge_id_type;

    /**
     * \brief The set of edges that are traversed during gather and scatter
     * operations.
     */
    enum edge_dir_type {
        /**
         * \brief No edges implies that no edges are processed during the
         * corresponding gather or scatter phase, essentially skipping
         * that phase.
         */
                NO_EDGES = 0,
        /**
         * \brief In edges implies that only whose target is the center
         * vertex are processed during gather or scatter.
         */
                IN_EDGES = 1,
        /**
         * \brief Out edges implies that only whose source is the center
         * vertex are processed during gather or scatter.
         */
                OUT_EDGES = 2 ,
        /**
         * \brief All edges implies that all adges adjacent to a the
         * center vertex are processed on gather or scatter.  Note that
         * some neighbors may be encountered twice if there is both an in
         * and out edge to that neighbor.
         */
                ALL_EDGES = 3};
} // end of namespace tux

#endif //TUX_GRAPH_BASIC_TYPES_HPP
