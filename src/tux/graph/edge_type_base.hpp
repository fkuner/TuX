//
// Created by fkuner on 2020/3/26.
//

#ifndef TUX_EDGE_TYPE_BASE_HPP
#define TUX_EDGE_TYPE_BASE_HPP

#include "graph_basic_types.hpp"

namespace tux {

    class edge_type_base {
    public:
        typedef tux::lvid_type lvid_type;
        typedef tux::edge_id_type edge_id_type;

        edge_type_base() {}

        edge_type_base(edge_id_type meid, lvid_type msource, lvid_type mtarget) {
            eid = meid;
            source = msource;
            target = mtarget;
        }

        lvid_type get_source(){ return source;}

        lvid_type get_target(){ return target;}

        edge_id_type get_eid(){ return eid; }
    private:
        edge_id_type eid;
        int partition_id = 0;
        lvid_type source;
        lvid_type target;

    };

}

#endif //TUX_EDGE_TYPE_BASE_HPP
