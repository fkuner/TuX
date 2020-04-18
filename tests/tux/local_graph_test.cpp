//
// Created by fkuner on 2020/4/8.
//

#include <iostream>

#include <tux/graph/local_graph.hpp>
#include <tux/graph/vertex_type_base.hpp>
#include <tux/graph/edge_type_base.hpp>
#include <tux/graph/graph_basic_types.hpp>


#include <graphlab/util/random.hpp>
#include <graphlab/macros_def.hpp>

using tux::vertex_type_base;
using tux::edge_type_base;
using tux::local_graph;

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

typedef local_graph<FVertexType, SVertexType, EdgeType> graph;
int main()
{

    typedef tux::lvid_type lvid_type;

    // test construct function
    graph g1;
    graph *g2 = new graph(10,10);

    std::vector<lvid_type> flvid_array{2,5,6};
    std::vector<lvid_type> slvid_array{0,1,3,4,7};


    // test add vertex
    for(auto id : flvid_array)
    {
        FVertexType fvertex(id,0);
        g1.add_fvertex(fvertex);
    }
    for(auto id : slvid_array)
    {
        SVertexType svertex(id, 0);
        g1.add_svertex(svertex);
    }

    // test add edge
    g1.add_edge(EdgeType(0,2,1,0));
    g1.add_edge(EdgeType(1,2,0,1));
    g1.add_edge(EdgeType(2,2,3,2));
    g1.add_edge(EdgeType(3,2,4,3));
    g1.add_edge(EdgeType(4,6,7,4));
    g1.add_edge(EdgeType(5,6,4,5));
    g1.add_edge(EdgeType(6,5,4,6));


    std::cout << "numbers of fvertex:" << g1.num_fvertices() << std::endl;
    std::cout << "numbers of svertex:" << g1.num_svertices() << std::endl;
    std::cout << "numbers of edge:" << g1.num_edges() << std::endl;

    // traverse all edges from the fvertex_array or svertex_array
    std::vector<FVertexType> fvertex_array = g1.get_fvertex_array();
    std::vector<EdgeType> edge_array = g1.get_edge_array();
    for(auto fvertex : fvertex_array){
        for(auto index : fvertex.edge_index_array) {
            std::cout<<edge_array[index].get_edge_data()<<" ";
        }
    }
    std::cout<<std::endl;


    std::cout << "test local_graph successfully" << std::endl;

}