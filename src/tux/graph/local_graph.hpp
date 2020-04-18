//
// Created by fkuner on 2020/3/20.
//

#ifndef TUX_LOCAL_GRAPH_H
#define TUX_LOCAL_GRAPH_H


#include <cmath>

#include <string>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <queue>
#include <algorithm>
#include <functional>
#include <fstream>
#include <iostream>

#include "graph_basic_types.hpp"

namespace tux{

    /// FVertexType is the first type of vertex data, and SVertexType is the second.
    template<typename FVertexType, typename SVertexType, typename EdgeType>
    class local_graph{

    public:
        // CONSTRUCTORS ============================================================>

        /** Create an empty local_graph. */
        local_graph() : finalized(false){ }

        /** Create a local_graph with fnverts fvertices and snverts svertices. */
        local_graph(size_t fnverts, size_t snverts):
            fvertex_array(fnverts),
            svertex_array(snverts),
            finalized(false) { }

        /** \brief Resets the local_graph state */
        void clear(){
            fvertex_array.clear();
            svertex_array.clear();
            edge_array.clear();
        }

        /** \brief Finalize the local_graph structure */
        void finalize();

        /** \brief Get the number of vertices */
        size_t num_vertices() const
        {
            return fvertex_array.size() + svertex_array.size();
        }

        /** \brief Get the number of fvertices */
        size_t num_fvertices() const{
            return fvertex_array.size();
        }

        /** \brief Get the number of svertices */
        size_t num_svertices() const
        {
            return svertex_array.size();
        }

        /** \brief Get the number of edges. */
        size_t num_edges() const
        {
            return edge_array.size();
        }

        /** Change the capacity of fvertices. */
        void reserve_farray(size_t num_fvertices)
        {
            ASSERT_GE(num_fvertices, fvertex_array.size());
            fvertex_array.reserve(num_fvertices);

        }

        /** Change the capacity of svertices */
        void reserve_sarray(size_t num_svertices)
        {
            ASSERT_GE(num_svertices, svertex_array.size());
            svertex_array.reserve(num_svertices);
        }

        /** \brief Add additional fvertices up to provided num_fvertices.
         * This will fail if resizing dowm.
         */
        void resize_farray(size_t num_fvertices)
        {
             ASSERT_GE(num_fvertices, fvertex_array.size());
             fvertex_array.resize(num_fvertices);
        }

        /** \brief Add additional svertices up to provided num_svertices.
         *  This will fail if resizing down.
         */
         void resize_sarray(size_t num_svertices)
        {
             ASSERT_GE(num_svertices, svertex_array.size());
             svertex_array.resize(num_svertices);
        }

        /** \brief Creates a vertex of the frsit type. Vertex ids are
         *  assigned in increasing order with the first vertex having
         *  id 0.
         */
        void add_fvertex(FVertexType& fvertex = FVertexType()){
//            lvid_type vid = fvertex.get_lvid();
//            if(vid >= fvertex_array.size()){
//                // Enable capacity doubling if resizing beyond capacity
//                if(vid >= fvertex_array.capacity()){
//                    const size_t new_size = std::max(2 * fvertex_array.capacity(), size_t(vid));
//                    fvertex_array.reserve(new_size);
//                }
//                fvertex_array.resize(vid+1);
//            }
            //fvertex_array[vid] = fvertex;
            fvertex_array.push_back(fvertex);
            flvid2index.insert(std::make_pair(fvertex.get_lvid(), fvertex_array.size()-1));
        }

        /** \brief Creates a vertex of the second type. Vertex ids are
         *  assigned in increasing order with the first vertex having
         *  id 0.
         */
        void add_svertex(SVertexType& svertex = SVertexType() ){
//            lvid_type vid = svertex.get_lvid();
//            if(vid >= svertex_array.size()){
//                // Enable capacity doubling if resizing beyond capacity
//                if(vid >= svertex_array.capacity()){
//                    const size_t new_size = std::max(2 * svertex_array.capacity(), size_t(vid));
//                    svertex_array.reserve(new_size);
//                }
//                svertex_array.resize(vid+1);
//            }
            //svertex_array[vid] = svertex;
            svertex_array.push_back(svertex);
            slvid2index.insert(std::make_pair(svertex.get_lvid(), fvertex_array.size()-1));
        }

        /** Create an edge connecting vertex source to vertex target. Any existing
         *  data will be cleared. Should not be called after finalization.
         */
        edge_id_type add_edge(EdgeType edge){
            if (finalized) {
                //logstream(LOG_FATAL)
                //        << "Attempting add edges to a finalized local_graph." << std::endl;
            }



            lvid_type source = edge.get_source();
            lvid_type target = edge.get_target();

            edge_id_type index = edge_array.size();

            if(flvid2index.find(source)!=flvid2index.end()){
                fvertex_array[flvid2index[source]].edge_index_array.push_back(index);
                svertex_array[flvid2index[target]].edge_index_array.push_back(index);
            } else{
                fvertex_array[flvid2index[target]].edge_index_array.push_back(index);
                svertex_array[flvid2index[source]].edge_index_array.push_back(index);
            }

            // Add the edge to the set fo edges
            edge_array.push_back(edge);

            // This is not the final edge_id, so we always return 0.
            return 0;
        }

        void add_edges(const std::vector<EdgeType> edges);

        /** \brief Return a fvertex of given ID */
        FVertexType fvertex(lvid_type vid)
        {
            //ASSERT_LT(vid, fvertex_array.size());
            lvid_type index = flvid2index[vid];
            return fvertex_array[index];
        }

        /** \brief Return a svertex of given ID */
        SVertexType svertex(lvid_type vid)
        {
            //ASSERT_LT(vid, svertex_array.size());
            lvid_type index = slvid2index[vid];
            return svertex_array[index];
        }

//        void display()
//        {
//            for(int i = 0; i < fvertex_array.end() - fvertex_array.begin();i++)
//                std::cout<<fvertex_array[i];
//            //
//        }

        const std::vector<FVertexType>& get_fvertex_array()
        {
            return fvertex_array;
        }

        const std::vector<SVertexType>& get_svertex_array()
        {
            return svertex_array;
        }

        const std::vector<EdgeType>& get_edge_array()
        {
            return edge_array;
        }

    private:
        std::vector<FVertexType> fvertex_array;
        std::vector<SVertexType> svertex_array;

        std::unordered_map<lvid_type, int> flvid2index;
        std::unordered_map<lvid_type, int> slvid2index;

        std::vector<EdgeType> edge_array;
        bool finalized;

        friend class local_graph_test;
    };

}



#endif //TUX_LOCAL_GRAPH_H