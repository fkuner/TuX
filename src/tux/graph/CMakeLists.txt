project(TuX)

add_executable(tmp tmp.cpp ingress/distributed_ingress_base.hpp
        ingress/distributed_hybrid_ingress.hpp
         builtin_parsers.hpp
        ingress/ingress_edge_decision.hpp )