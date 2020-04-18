//
// Created by fkuner on 2020/4/9.
//

#include <iostream>
#include <tux.hpp>
#include <graphlab/rpc/dc.hpp>
#include <graphlab/util/mpi_tools.hpp>
#include <graphlab/rpc/dc_init_from_mpi.hpp>
#include <tux/graph/distributed_graph.hpp>
// #include <google/malloc_extension.h>
#include <graphlab/macros_def.hpp>
#include <tux/graph/vertex_type_base.hpp>
#include <tux/graph/edge_type_base.hpp>
#include <tux/graph/graph_basic_types.hpp>

using tux::vertex_type_base;
using tux::edge_type_base;

class FVertexType:public vertex_type_base{
private:

    int fvertex_data;
public:

    FVertexType(): vertex_type_base(){ }

    FVertexType(int _fvertex_data):vertex_type_base(){
        fvertex_data = _fvertex_data;
    }

    FVertexType(lvid_type _lvid, int _fvertex_data):vertex_type_base(_lvid){
        fvertex_data = _fvertex_data;
    }
};

class SVertexType:public vertex_type_base{
private:

    int svertex_data;
public:

    SVertexType(): vertex_type_base(){ }

    SVertexType(int _svertex_data):vertex_type_base(){
        svertex_data = _svertex_data;
    }

    SVertexType(lvid_type _lvid, int _svertex_data):vertex_type_base(_lvid){
        svertex_data = _svertex_data;
    }
};

class EdgeType:public edge_type_base{
private:
    int edge_data;
public:

    EdgeType(){

    }

    EdgeType(int _edge_data)
    {
        edge_data = _edge_data;
    }

    EdgeType(edge_id_type _eid, lvid_type _source, lvid_type _target, int _edge_data):
            edge_type_base(_eid,_source, _target){
        edge_data = _edge_data;
    }

    int get_edge_data()
    {
        return edge_data;
    }
};

typedef tux::distributed_graph<FVertexType, SVertexType, EdgeType> graph_type;

int main(int argc, char** argv)
{
    /// initial
}