//
// Created by fkuner on 2020/3/28.
//

#ifndef TUX_DISTRIBUTED_GRAPH_HPP
#define TUX_DISTRIBUTED_GRAPH_HPP

#ifndef __NO_OPENMP__
#include <omp.h>
#endif

#include <cmath>

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <graphlab/util/dense_bitset.hpp>


#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>
#include <sstream>


#include <tux/graph/ingress/distributed_ingress_base.hpp>
#include <tux/graph/ingress/distributed_hybrid_ingress.hpp>

#include <tux/graph/local_graph.hpp>


#include <graphlab/graph/graph_hash.hpp>

#include <graphlab/rpc/dc.hpp>
#include <graphlab/rpc/dc_dist_object.hpp>
#include <graphlab/rpc/buffered_exchange.hpp>

#include <graphlab/options/graphlab_options.hpp>


#include <graphlab/util/hopscotch_map.hpp>
#include <graphlab/graph/vertex_set.hpp>

using graphlab::distributed_control;
using graphlab::graphlab_options;
using graphlab::dc_dist_object;
using graphlab::procid_t;
using graphlab::hopscotch_map;
using graphlab::vertex_set;
using graphlab::iarchive;
using graphlab::oarchive;
using namespace graphlab::graph_hash;

#include <graphlab/macros_def.hpp>

namespace tux{
    template<typename FVertexType, typename SVertexType, typename EdgeType>
    class distributed_graph{
    public:

        typedef boost::function<bool(distributed_graph&, const std::string&, const std::string&)> line_parser_type;

        typedef graphlab::fixed_dense_bitset<RPC_MAX_N_PROCS> mirror_type;

        /// The type of the local graph used to store the graph data
//#ifdef USE_DYNAMIC_LOCAL_GRAPH
//        typedef tux::dynamic_local_graph<VertexData, EdgeData> local_graph_type;
//#else
        typedef tux::local_graph<FVertexType, SVertexType, EdgeType> local_graph_type;
//#endif
        typedef tux::distributed_graph<FVertexType, SVertexType, EdgeType> graph_type;

        //friend class distributed_ingress_base<FVertexType, SVertexType, EdgeType>;

        // Make friends with Ingress classes
        //friend class distributed_hybrid_ingress<FVertexType, SVertexType, EdgeType>

        typedef tux::vertex_id_type vertex_id_type;
        typedef tux::lvid_type lvid_type;
        typedef tux::edge_id_type edge_id_type;


        // CONSTRUCTORS ==========================================================>
        distributed_graph(distributed_control& dc, const graphlab_options& opts = graphlab_options()):
            rpc(dc, this), finalized(false), vid2lvid(),
            nverts(0), nedges(0), local_own_nverts(0), nreplicas(0),
            ingress_ptr(NULL),
#ifdef _OPENMP
            vertex_exchange(dc, omp_get_max_threads()),
#else
            vertex_exchange(dc),
#endif
            vset_exchange(dc), parallel_ingress(true){
            rpc.barrier();
            set_options(opts);
        };

        ~distributed_graph(){
            delete ingress_ptr;
            ingress_ptr = NULL;
        }

    private:
        void set_options(const graphlab_options& opts)
        {
            std::cout<<"not implemented";
        }
    public:

        // METHODS ===============================================================>
        void finalize(){

        }

        bool is_finalized(){
            return finalized;
        }

        /** \brief Get the number of vertices */
        size_t num_vertices() const {return nverts;}

        /** \brief Get the number of edges */
        size_t num_edges() const { return nedges; }

        FVertexType fvertex(vertex_id_type vid) {
            return vertex_type(*this, local_vid(vid));
        }

        size_t num_in_edges(const vertex_id_type vid) const {
            return get_vertex_record(vid).num_in_edges;
        }

        size_t num_out_edges(const vertex_id_type vid) const {
            return get_vertex_record(vid).num_out_edges;
        }

        bool add_fvertex(const vertex_id_type& vid,
                        const FVertexType& fvdata = FVertexType() ) {
        }

        bool add_svertex(const vertex_id_type& vid,
                         const SVertexType& svdata = SVertexType() ) {
        }

        bool add_edge(vertex_id_type source, vertex_id_type target,
                      const EdgeType& edata = EdgeType()) {
        }

        template <typename TransformType>
        void transform_vertices(TransformType transform_functor,
                                const vertex_set vset = complete_set()) {

        }

        template <typename TransformType>
        void transform_edges(TransformType transform_functor,
                             const vertex_set& vset = complete_set(),
                             edge_dir_type edir = IN_EDGES) {

        }

        /** \brief Load the graph from an archive */
        void load(iarchive& arc) {
            // read the vertices
            arc >> nverts
                >> nedges
                >> local_own_nverts
                >> nreplicas
                >> vid2lvid
                >> lvid2record
                >> local_graph;
            finalized = true;
            // check the graph condition
        } // end of load

        /** \brief Save the graph to an archive */
        void save(oarchive& arc) const {
            if(!finalized) {
                logstream(LOG_FATAL)
                    << "\n\tAttempting to save a graph before calling graph.finalize()."
                    << std::endl;
            }
            // Write the number of edges and vertices
            arc << nverts
                << nedges
                << local_own_nverts
                << nreplicas
                << vid2lvid
                << lvid2record
                << local_graph;
        } // end of save

        /// \brief Clears and resets the graph, releasing all memory used.
        void clear () {
        }

        bool load_binary(const std::string& prefix) {
        };

        bool save_binary(const std::string& prefix) {

        }

        template<typename Writer>
        void save_to_posixfs(const std::string& prefix, Writer writer,
                             bool gzip = true,
                             bool save_vertex = true,
                             bool save_edge = true,
                             size_t files_per_machine = 4) {
        }

        template<typename Writer>
        void save_to_hdfs(const std::string& prefix, Writer writer,
                          bool gzip = true,
                          bool save_vertex = true,
                          bool save_edge = true,
                          size_t files_per_machine = 4) {

        }

        template<typename Writer>
        void save(const std::string& prefix, Writer writer,
                  bool gzip = true, bool save_vertex = true, bool save_edge = true,
                  size_t files_per_machine = 4) {

        }

        void save_format(const std::string& prefix, const std::string& format,
                         bool gzip = true, size_t files_per_machine = 4) {

        }

        void load_from_posixfs(std::string prefix,
                               line_parser_type line_parser) {

        }

        void load_from_hdfs(std::string prefix, line_parser_type line_parser) {

        }

        void load(std::string prefix, line_parser_type line_parser) {

        } // end of load

        void load_synthetic_powerlaw(size_t nverts, bool in_degree = false,
                                     double alpha = 2.1, size_t truncate = (size_t)(-1)) {
        }

        void load_format(const std::string& path, const std::string& format) {

        }

        /****************************************************************************
         *                     Vertex Set Functions                                 *
         *                     ----------------------                               *
         * Manages operations involving sets of vertices                            *
         ****************************************************************************/

        /**
          * \brief Retuns an empty set of vertices
          */
        static vertex_set empty_set() {
            return vertex_set(false);
        }

        /**
         *  \brief Retuns a full set of vertices
         */
        static vertex_set complete_set() {
            return vertex_set(true);
        }


        void sync_vertex_set_master_to_mirrors(vertex_set& vset) {

        }

        size_t vertex_set_size(const vertex_set& vset) {

        }

        bool vertex_set_empty(const vertex_set& vset) {

        }


        /**
         * \internal
         * The vertex record stores information associated with each
         * vertex on this proc
         */
        struct vertex_record {
            /// The official owning processor for this vertex
            procid_t owner;
            /// The local vid of this vertex on this proc
            vertex_id_type gvid;
            /// The number of in edges
            vertex_id_type num_in_edges, num_out_edges;
            /** The set of proc that mirror this vertex.  The owner should
                NOT be in this set.*/
            mirror_type _mirrors;
            vertex_record() :
                    owner(-1), gvid(-1), num_in_edges(0), num_out_edges(0) { }
            vertex_record(const vertex_id_type& vid) :
                    owner(-1), gvid(vid), num_in_edges(0), num_out_edges(0) { }
            procid_t get_owner () const { return owner; }
            const mirror_type& mirrors() const { return _mirrors; }
            size_t num_mirrors() const { return _mirrors.popcount(); }

            void clear() {
                _mirrors.clear();
            }

            void load(iarchive& arc) {
                clear();
                arc >> owner
                    >> gvid
                    >> num_in_edges
                    >> num_out_edges
                    >> _mirrors;
            }

            void save(oarchive& arc) const {
                arc << owner
                    << gvid
                    << num_in_edges
                    << num_out_edges
                    << _mirrors;
            } // end of save

            bool operator==(const vertex_record& other) const {
                return (
                        (owner == other.owner) &&
                        (gvid == other.gvid)  &&
                        (num_in_edges == other.num_in_edges) &&
                        (num_out_edges == other.num_out_edges) &&
                        (_mirrors == other._mirrors)
                );
            }
        }; // end of vertex_record

        /** \internal
         * \brief converts a local vertex ID to a local vertex object
         */
        local_vertex_type l_vertex(lvid_type vid) {
            return local_vertex_type(*this, vid);
        }

        /** \internal
     *\brief Get the Total number of vertex replicas in the graph */
        size_t num_replicas() const { return nreplicas; }

        /** \internal
         *\brief Get the number of vertices local to this proc */
        size_t num_local_vertices() const { return local_graph.num_vertices(); }

        /** \internal
         *\brief Get the number of edges local to this proc */
        size_t num_local_edges() const { return local_graph.num_edges(); }

        /** \internal
         *\brief Get the number of vertices owned by this proc */
        size_t num_local_own_vertices() const { return local_own_nverts; }

        /** \internal
         *\brief Convert a global vid to a local vid */
        lvid_type local_vid (const vertex_id_type vid) const {
            // typename boost::unordered_map<vertex_id_type, lvid_type>::
            //   const_iterator iter = vid2lvid.find(vid);
            typename hopscotch_map_type::const_iterator iter = vid2lvid.find(vid);
            return iter->second;
        } // end of local_vertex_id

        /** \internal
         *\brief Convert a local vid to a global vid */
        vertex_id_type global_vid(const lvid_type lvid) const {
            ASSERT_LT(lvid, lvid2record.size());
            return lvid2record[lvid].gvid;
        } // end of global_vertex_id

        /** \internal
     * \brief Returns true if the local graph as an instance of (master or mirror)
     * of the vertex ID.
     */
        bool contains_vertex(const vertex_id_type vid) const {
            return vid2lvid.find(vid) != vid2lvid.end();
        }
        /**
         * \internal
         * \brief Returns an edge list of all in edges of a local vertex ID
         *        on the local graph
         *
         * Equivalent to l_vertex(lvid).in_edges()
         */
        local_edge_list_type l_in_edges(const lvid_type lvid) {
            return local_edge_list_type(*this, local_graph.in_edges(lvid));
        }

        /**
         * \internal
         * \brief Returns the number of in edges of a local vertex ID
         *        on the local graph
         *
         * Equivalent to l_vertex(lvid).num in_edges()
         */
        size_t l_num_in_edges(const lvid_type lvid) const {
            return local_graph.num_in_edges(lvid);
        }

        /**
         * \internal
         * \brief Returns an edge list of all out edges of a local vertex ID
         *        on the local graph
         *
         * Equivalent to l_vertex(lvid).out_edges()
         */
        local_edge_list_type l_out_edges(const lvid_type lvid) {
            return local_edge_list_type(*this, local_graph.out_edges(lvid));
        }

        /**
         * \internal
         * \brief Returns the number of out edges of a local vertex ID
         *        on the local graph
         *
         * Equivalent to l_vertex(lvid).num out_edges()
         */
        size_t l_num_out_edges(const lvid_type lvid) const {
            return local_graph.num_out_edges(lvid);
        }

        procid_t procid() const {
            return rpc.procid();
        }


        procid_t numprocs() const {
            return rpc.numprocs();
        }

        distributed_control& dc() {
            return rpc.dc();
        }

        const vertex_record& get_vertex_record(vertex_id_type vid) const {
            // typename boost::unordered_map<vertex_id_type, lvid_type>::
            //   const_iterator iter = vid2lvid.find(vid);
            typename hopscotch_map_type::const_iterator iter = vid2lvid.find(vid);
            ASSERT_TRUE(iter != vid2lvid.end());
            return lvid2record[iter->second];
        }

        /** \internal
         * \brief Returns the internal vertex record of a given local vertex ID
         */
        vertex_record& l_get_vertex_record(lvid_type lvid) {
            ASSERT_LT(lvid, lvid2record.size());
            return lvid2record[lvid];
        }

        /** \internal
         * \brief Returns the internal vertex record of a given local vertex ID
         */
        const vertex_record& l_get_vertex_record(lvid_type lvid) const {
            ASSERT_LT(lvid, lvid2record.size());
            return lvid2record[lvid];
        }

        /** \internal
         * \brief Returns true if the provided global vertex ID is a
         *        master vertex on this machine and false otherwise.
         */
        bool is_master(vertex_id_type vid) const {
            const procid_t owning_proc = hash_vertex(vid) % rpc.numprocs();
            return (owning_proc == rpc.procid());
        }


        procid_t master(vertex_id_type vid) const {
            const procid_t owning_proc = hash_vertex(vid) % rpc.numprocs();
            return owning_proc;
        }

        /** \internal
         * \brief Returns true if the provided local vertex ID is a master vertex.
         *        Returns false otherwise.
         */
        bool l_is_master(lvid_type lvid) const {
            ASSERT_LT(lvid, lvid2record.size());
            return lvid2record[lvid].owner == rpc.procid();
        }

        /** \internal
         * \brief Returns the master procid for vertex lvid.
         */
        procid_t l_master(lvid_type lvid) const {
            ASSERT_LT(lvid, lvid2record.size());
            return lvid2record[lvid].owner;
        }


        /** \internal
         *  \brief Returns a reference to the internal graph representation
         */
        local_graph_type& get_local_graph() {
            return local_graph;
        }

        /** \internal
         *  \brief Returns a const reference to the internal graph representation
         */
        const local_graph_type& get_local_graph() const {
            return local_graph;
        }




        /** \internal
         * This function synchronizes the master vertex data with all the mirrors.
         * This function must be called simultaneously by all machines
         */
        void synchronize(const vertex_set& vset = complete_set()) {
            typedef std::pair<vertex_id_type, vertex_data_type> pair_type;

            procid_t sending_proc;
            // Loop over all the local vertex records


#ifdef _OPENMP
#pragma omp parallel for
#endif
            for(lvid_type lvid = 0; lvid < lvid2record.size(); ++lvid) {
                typename buffered_exchange<pair_type>::buffer_type recv_buffer;
                const vertex_record& record = lvid2record[lvid];
                // if this machine is the owner of a record then send the
                // vertex data to all mirrors
                if(record.owner == rpc.procid() && vset.l_contains(lvid)) {
                    foreach(size_t proc, record.mirrors()) {
                        const pair_type pair(record.gvid, local_graph.vertex_data(lvid));
#ifdef _OPENMP
                        vertex_exchange.send(proc, pair, omp_get_thread_num());
#else
                        vertex_exchange.send(proc, pair);
#endif
                    }
                }
                // Receive any vertex data and update local mirrors
                while(vertex_exchange.recv(sending_proc, recv_buffer, true)) {
                    foreach(const pair_type& pair, recv_buffer)  {
                        vertex(pair.first).data() = pair.second;
                    }
                    recv_buffer.clear();
                }
            }

            typename buffered_exchange<pair_type>::buffer_type recv_buffer;
            vertex_exchange.flush();
            while(vertex_exchange.recv(sending_proc, recv_buffer)) {
                foreach(const pair_type& pair, recv_buffer) {
                    vertex(pair.first).data() = pair.second;
                }
                recv_buffer.clear();
            }
            ASSERT_TRUE(vertex_exchange.empty());
        } // end of synchronize

    private:
        mutable dc_dist_object<distributed_graph> rpc;

        bool finalized;

        /** The local graph data */
        local_graph_type local_graph;

        /** The map from global vertex ids to vertex records */
        std::vector<vertex_record>  lvid2record;

        // boost::unordered_map<vertex_id_type, lvid_type> vid2lvid;
        /** The map from global vertex ids back to local vertex ids */
        typedef hopscotch_map<vertex_id_type, lvid_type> hopscotch_map_type;
        typedef hopscotch_map_type vid2lvid_map_type;

        hopscotch_map_type vid2lvid;

        /** The global number of vertices and edges */
        size_t nverts, nedges;

        /** The number of vertices owned by this proc */
        size_t local_own_nverts;

        /** The global number of vertex replica */
        size_t nreplicas;

        /** pointer to the distributed ingress object*/
        distributed_ingress_base<FVertexType, SVertexType, EdgeType>* ingress_ptr;

        /** Buffered Exchange used by synchronize() */
        buffered_exchange<std::pair<vertex_id_type, vertex_data_type> > vertex_exchange;

        /** Buffered Exchange used by vertex sets */
        buffered_exchange<vertex_id_type> vset_exchange;

        /** Command option to disable parallel ingress. Used for simulating single node ingress */
        bool parallel_ingress;

        lock_manager_type lock_manager;

        void set_ingress_method(const std::string& method,
                                size_t bufsize = 50000, bool usehash = false, bool userecent = false) {

        }
    };
}

#include <graphlab/macros_undef.hpp>

#endif //TUX_DISTRIBUTED_GRAPH_HPP
