//
// Created by fkuner on 2020/3/26.
//

#ifndef TUX_VERTEX_TYPE_BASE_HPP
#define TUX_VERTEX_TYPE_BASE_HPP

#include <vector>
#include "graph_basic_types.hpp"

namespace tux{

    class vertex_type_base{
    public:
        typedef tux::lvid_type lvid_type;
        typedef tux::edge_id_type edge_ie_type;

        vertex_type_base(){ }

        vertex_type_base(lvid_type mlvid){
            lvid = mlvid;
        }

        vertex_type_base(lvid_type mlvid, const std::vector<edge_id_type> &medge_index_array)
        {
            lvid = mlvid;
            edge_index_array.assign(medge_index_array.begin(), medge_index_array.end());
        }

        std::vector<edge_id_type> edge_index_array;

        lvid_type get_lvid()
        {
            return lvid;
        }

    private:

        /// \brief Return the ID of the vertex.
        lvid_type lvid;


    };
}

#endif //TUX_VERTEX_TYPE_BASE_HPP
