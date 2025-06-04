#include "graph.h"
#include "reduction.h"
#include "allocation.h"
#include "extraction.h"
#include "report.h"
#include "verification.h"
#include "flownetwork.h"
#include "lp.h"
#include <iostream>
#include "args.h"
#include "ratioselection.h"
#include "app.h"
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <set>

using Heap = boost::heap::fibonacci_heap<std::pair<int, VertexID>>;
int main(int argc, char **argv) {
    Args args = Args();
    args.argsParse(argc, argv);
    bool is_directed = args.getOption("-t") == "d";
    bool is_exact = args.getOption("-a") == "e";
    bool is_vw = args.getOption("-vw") == "t";
    bool is_exp = args.getOption("-exp") == "t";
    bool is_dc = args.getOption("-dc") == "t";
    bool is_seq = args.getOption("-seq") == "t";
    bool is_reduction_ablation = args.getOption("-ra") == "t";
    bool is_map = args.getOption("-map") == "t";
    bool is_res = args.getOption("-res") == "t";
    bool is_mul = args.getOption("-multi") == "t";
    bool stable = args.getOption("-stable") == "t";
    bool e_stable = args.getOption("-estable") == "t";
    bool is_stats = args.getOption("-stats") == "t";
    bool is_sample = args.getOption("-sample") == "t";
    bool is_print_c = args.getOption("-printc") == "t";
    bool is_wcore_shrink = args.getOption("-wshrink") == "t";
    bool is_init_wcore = args.getOption("-initwcore") == "t";
    bool is_original = args.getOption("-original") == "t";
    bool is_log = args.getOption("-log") == "t";

    double sample_rate = std::stod(args.getOption("-rate"));
    double epsilon = std::stod(args.getOption("-eps"));
    double learning_rate = std::stod(args.getOption("-gamma"));

    ui order_type = std::stoi(args.getOption("-o"));
    ui iter_num = std::stoi(args.getOption("-it"));
    ui res_width = std::stoi(args.getOption("-width"));
//    int parallel_thread_num = std::stoi(args.getOption("-p"));
//    std::string ompNumThreads = "OMP_NUM_THREADS=";
//    ompNumThreads += std::to_string(parallel_thread_num);
//
//    char* env = new char[ompNumThreads.size() + 1];
//    std::strcpy(env, ompNumThreads.c_str());
//    putenv(env);

    std::string red_type = args.getOption("-red");
    std::string alloc_type = args.getOption("-alloc");
    std::string ext_type = args.getOption("-ext");
    std::string ver_type = args.getOption("-ver");
    //todo
    Graph graph = Graph(is_directed);
    if(!is_sample) graph.loadGraphFromFile(args.getOption("-path"));
    else graph.loadGraphFromFile(args.getOption("-path"), is_sample, sample_rate);
    graph.removeMultiEdges(graph);
//    graph.subgraph_density_upper_bound = 1e20;
    printf("io finished\nn:%u m:%u\n", graph.getVerticesCount(), graph.getEdgesCount());
    auto begin = std::chrono::steady_clock::now();
    graph.init();
    if(is_sample) graph.sample(sample_rate);
    Reduction red;
    Allocation alloc;
    Extraction ext;
    Verification ver;
    Report rep;
    rep.add_graph_size(graph.getVerticesCount(), graph.getEdgesCount());
    if(iter_num == -1){
        if(graph.getEdgesCount()<1e7) iter_num = 100;
        else iter_num = 10;
        printf("Default iter_num: %u\n", iter_num);
    }
    if(!is_exact&&!is_original&&graph.getEdgesCount()<1e7){
        if(is_init_wcore){
            printf("The graph is too small for w-core init in appro, set is_init_wcore = false\n");
            is_init_wcore = false;
        }
    }
    auto wcore_init_density = [&rep](Graph &graph){
        printf("Start w-core initialization\n");
        auto begin_init_wcore = std::chrono::steady_clock::now();
        WCore w_core_;
        Graph tmp(true, 0);
        w_core_.generateInitWCore(graph,tmp);
        graph.subgraph_density=std::max(tmp.subgraph_density,graph.subgraph_density);
        printf("Begin with w-core density %.10lf\n", graph.subgraph_density);
        auto end_init_wcore = std::chrono::steady_clock::now();
        printf("Init w-core time: %f\n", std::chrono::duration<double>(end_init_wcore - begin_init_wcore).count());
        rep.add_total_wcore_time(std::chrono::duration<double>(end_init_wcore - begin_init_wcore).count());
    };
    if (!is_exact) {
        if (!is_directed) {
            if(red_type == "k-core" || red_type == "stable") graph.coreReduce(graph);
            double ratio;
            bool flag = true;
            bool lp_type = alloc_type == "fista"? 1: 0;
            LinearProgramming lp(is_directed, lp_type, 0, 0, order_type);
            if(alloc_type == "fw" || alloc_type == "fista" || alloc_type == "mwu"){
                if(order_type == 0){

                }
                if(order_type == 1) graph.coreOrder(graph);
                if(order_type == 2) graph.CoreOrder(graph);
                lp.Init(graph);
            }
            CoreApp ca = CoreApp();
            PKMC pkmc = PKMC();
            Greedy gr = Greedy();
            if(alloc_type == "greedypp") gr.Init(graph);
            if(alloc_type == "core-app"){
                ca.Init(graph);
            }
            if(alloc_type == "pkmc") pkmc.Init(graph);
            auto vertices = new std::vector<VertexID>[1];
            ui T = 1;
            auto end = std::chrono::steady_clock::now();
            double opt = 0;
            graph.subgraph_density_upper_bound = 1e9;
            while (flag) {
                std::cout<<"T = "<<T<<std::endl;
                if(alloc_type == "flow-app" && red_type == "k-core") graph.coreReduce(graph, (ui) (graph.subgraph_density_lower_bound));
                if(red_type == "k-core" && (alloc_type == "fw" || alloc_type == "fista" || alloc_type == "mwu")) red.UndirectedkCoreReduction(graph, lp);
                if(red_type == "stable" && (alloc_type == "fw" || alloc_type == "fista" || alloc_type == "mwu")) red.UndirectedStableReduction(graph, lp);
                if(alloc_type == "core-app") alloc.UndirectedCoreAppAllocation(graph, ca);
                if(alloc_type == "greedy"){
                    alloc.UndirectedGreedyAllocation(graph);
                    break;
                }
                if(alloc_type == "flow-app"){
                    alloc.UndirectedFlowAppAllocation(graph);
                }
                if(alloc_type == "greedypp") alloc.UndirectedGreedyppAllocation(graph, gr);
                if(alloc_type == "fw") alloc.UndirectedlpAllocation(graph, lp, T, is_seq); 
                if(alloc_type == "fista") alloc.UndirectedFistaAllocation(graph, lp, T, is_seq); 
                if(alloc_type == "mwu") alloc.UndirectedMWUAllocation(graph, lp, T, is_seq);
                if(ext_type == "cp") ext.UndirectedlpExactExtraction(graph, lp, vertices);
                if(ver_type == "cp") flag = ver.UndirectedLpAppVerification(graph, lp, vertices, epsilon);
                if(ver_type == "core-app") flag = ver.UndirectedCoreAppVerification(graph, ca);
                if(ver_type == "flow-app") flag = ver.UndirectedFlowAppVerification(graph, epsilon);
                T<<=1;
                end = std::chrono::steady_clock::now();
                if(graph.subgraph_density_lower_bound > opt){
                    opt = graph.subgraph_density_lower_bound;
                    printf("%.10lf %.10lf\n", opt, std::chrono::duration<double>(end - begin).count());
                }
            }
        } else if(!is_vw){
            std::pair<double, double> ratio;
            double ratio_o = 0, ratio_p = 0;
            double reduction_ratio = 0;
            double total_vertices_num = 0;
            double total_edges_num = 0;
            double total_iter_num = 0;
            RatioSelection ratio_selection(graph);
            bool is_init_ratio = false;
            ui ratio_count = 0;
            if(is_init_wcore){
                wcore_init_density(graph);
            }
            while (ratio_selection.ratioSelection(graph.getVerticesCount(),
                                                  ratio,
                                                  is_init_ratio,
                                                  is_vw,
                                                  is_dc,
                                                  ratio_o,
                                                  ratio_p,
                                                  graph.subgraph_density,
                                                  epsilon, is_res, res_width, is_map)) {
                ratio_count++;
                ui edges_count = 0;
                bool flag = true;
                bool is_init_red = false;
                bool is_init_lp = false;
                bool is_reduced = false;
                bool is_stable_set = false;
                double rho, rho_c;
                double l = learning_rate * graph.subgraph_density;
                double r = graph.subgraph_density_upper_bound;
                FlowNetwork flow;
                WCore w_core;
                int lp_type= alloc_type == "fista"? 1: 0;
                //LinearProgramming lp(is_directed, 0, 0, 0, order_type);
                LinearProgramming lp(is_directed, lp_type, 0, 0, order_type);
                std::vector<std::pair<VertexID, VertexID>> edges;
                std::vector<std::vector<VertexID>> vertices(2);
                std::vector<Heap> heap;
                std::vector<std::vector<Heap::handle_type>> handles;
                std::vector<std::vector<bool>> is_peeled;
                std::pair<ui, ui> best_pos(0, 0);
                static Graph subgraph(is_directed, 0);
                subgraph.subgraph_density = graph.subgraph_density;
                rep.add_density(total_iter_num, graph.subgraph_density);

                if (red_type != "appro-xy-core" && red_type != "exact-xy-core"){
                    subgraph = Graph(is_directed,0); // For core inheritance
                }

                if (is_wcore_shrink) {
                    static double last_density = 0;
                    if (graph.subgraph_density > last_density * 1.25) {
                        last_density = graph.subgraph_density;
                        
                        auto begin_wcore = std::chrono::steady_clock::now();
                        if(is_log)
                            printf("Enter Wcore, edges #: %d, density %.10lf\n", graph.getEdgesCount(), graph.subgraph_density);

                        Graph _wcore = Graph(is_directed, 0);
                        w_core.getWCore(graph, _wcore, ceil(graph.subgraph_density*graph.subgraph_density/4.0));
                        graph.deg_[0] = _wcore.deg_[0];
                        graph.deg_[1] = _wcore.deg_[1];
                        graph.adj_[0] = _wcore.adj_[0];
                        graph.adj_[1] = _wcore.adj_[1];
                        graph.edges_count_ = _wcore.edges_count_;
                        is_reduced = is_init_red = false;
                        
                        auto end_wcore = std::chrono::steady_clock::now();
                        if(is_log)
                            printf("Exit Wcore, edges #: %d\n", graph.getEdgesCount());

                        rep.add_graph_size(graph.getVerticesCount(), graph.getEdgesCount()); // Report
                        rep.add_total_wcore_time(std::chrono::duration<double>(end_wcore - begin_wcore).count());
                    }
                }

                while (flag) {
                    if (!is_reduced || (alloc_type != "fw" && alloc_type != "fista")) {
                        is_reduced = true;
                        auto begin_core = std::chrono::steady_clock::now();
                        if (red_type == "exact-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, is_map, true, true, is_res, res_width, false);
                        } else if (red_type == "appro-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, is_map, false, true, is_res, res_width, false);
                        } else if (red_type == "w-core") {
                            red.wCoreReduction(graph, subgraph, w_core);
                        } else {
                            subgraph = graph;
                        }
                        if (is_stats) {
                            total_vertices_num += subgraph.getVerticesCount();
                            total_edges_num += subgraph.getEdgesCount();
                        }
                        if (subgraph.getEdgesCount() == 0) {
                            double c;
                            if (is_map) {
                                if (ratio.first < 1 && ratio.second > 1) {
                                    c = 1;
                                } else if (ratio.second <= 1) {
                                    c = (ratio.first + ratio.second) / 2;
                                } else if (ratio.first >= 1) {
                                    c = 2 / (1 / ratio.first + 1 / ratio.second);
                                }
                            } else
                                c = (ratio.first + ratio.second) / 2;
                            if (is_res) {
                                ratio_o = std::max(ratio.first, c / res_width);
                                ratio_p = std::min(ratio.second, c * res_width);
                            } else {
                                ratio_o = ratio.first;
                                ratio_p = ratio.second;
                            }
                            break;
                        }
                        auto end_core = std::chrono::steady_clock::now();
                        rep.add_total_xycore_time(std::chrono::duration<double>(end_core - begin_core).count());
                    }
                    if (is_mul && (alloc_type == "fw"||alloc_type=="fista") && subgraph.subgraph_density < graph.subgraph_density) {
                        auto begin_red = std::chrono::steady_clock::now();
                        subgraph.subgraph_density = graph.subgraph_density;
                        is_init_red = false;
                        red.xyCoreReduction(subgraph, subgraph, ratio, subgraph.subgraph_density, r, is_init_red, is_dc, is_map, false, false, is_res, res_width, true);
                        auto end_red = std::chrono::steady_clock::now();
                        rep.add_total_xycore_time(std::chrono::duration<double>(end_red - begin_red).count());
                    }
                    if (alloc_type == "fw"||alloc_type=="fista") 
                        red.stableSetReduction(subgraph, lp, edges, is_stable_set, true);
                    auto begin_iter = std::chrono::steady_clock::now();
                    if (alloc_type == "greedy")
                        alloc.directedBSApproAllocation(graph, ratio, heap, handles, is_peeled, edges_count, is_init_lp);
                    if (alloc_type == "xy-core-appro")
                        alloc.xyCoreApproAllocation(subgraph, best_pos);
                    if (alloc_type == "w-core-appro")
                        alloc.wCoreApproAllocation(subgraph, w_core, best_pos);
                    if (alloc_type == "fw")
                        alloc.directedCPAllocation(subgraph, lp, iter_num, is_init_lp, ratio, !is_seq, is_exp, is_map);
                    if (alloc_type == "fista")
                    {
                        alloc.directedFistaAllocation(subgraph, lp, iter_num, is_init_lp, ratio, !is_seq, is_exp, is_map, is_log);
                    }
                    auto end_iter = std::chrono::steady_clock::now();
                    rep.add_total_iter_time(std::chrono::duration<double>(end_iter - begin_iter).count());

                    if (ext_type == "core-appro")
                        ext.directedCoreApproExtraction(graph, subgraph, best_pos);
                    if (ext_type == "cp")
                    {
                        ext.directedCPExtraction(subgraph, lp, best_pos, vertices, ratio, ratio_o, ratio_p, rho, rho_c,
                                                 is_map);
                    }
                    if (ext_type == "frac")
                        ext.directedFractionalPeelingExtraction(subgraph, lp, best_pos, vertices, ratio, ratio_o, ratio_p, rho, rho_c,
                                                                is_map);
                    if (ext_type == "greedy")
                        ext.directedBSApproExtraction(graph, is_peeled, vertices);

                    auto begin_ver = std::chrono::steady_clock::now();
                    if (ver_type == "greedy")
                        flag = ver.directedBSApproVerification(graph, edges_count, vertices);
                    if (ver_type == "no")
                        break;
                    if (ver_type == "cp")
                        flag = ver.directedCPVerification(graph, subgraph, lp, best_pos, vertices, ratio, rho, rho_c,
                                                          ratio_o, ratio_p, is_stable_set, edges, epsilon, stable,
                                                          is_map);
                    auto end_ver = std::chrono::steady_clock::now();
                    rep.add_total_verify_time(std::chrono::duration<double>(end_ver - begin_ver).count());
                }
                if (ext_type == "core-appro")
                    break;

                total_iter_num += lp.cur_iter_num;
                rep.add_density(total_iter_num, graph.subgraph_density);

                double c;
                if (is_map) {
                    if (ratio.first < 1 && ratio.second > 1) {
                        c = 1;
                    } else if (ratio.second <= 1) {
                        c = (ratio.first + ratio.second) / 2;
                    } else if (ratio.first >= 1) {
                        c = 2 / (1 / ratio.first + 1 / ratio.second);
                    }
                } else
                    c = (ratio.first + ratio.second) / 2;
            }
            printf("ratio count: %d, density: %f, S/T: %lu/%lu\n", ratio_count, graph.subgraph_density,
                   graph.vertices[0].size(), graph.vertices[1].size());
            if (is_stats){
                printf ("avg vertices #: %f\navg edges #: %f\navg iterations #: %f\n", total_vertices_num / ratio_count, total_edges_num / ratio_count, total_iter_num / ratio_count);
                printf ("total subgraph percentage:\nvertex #: %.4lf edge #: %.4lf\n",
                            total_vertices_num * 1.0 / ratio_count / graph.getVerticesCount(), 
                            total_edges_num * 1.0 / ratio_count / graph.getEdgesCount());
            }
            rep.set_final_density(graph.subgraph_density);
            rep.set_total_iteration(total_iter_num);
        } else {
            std::pair<double, double> ratio;
            double ratio_o = 0, ratio_p = 0;
            double reduction_ratio = 0;
            RatioSelection ratio_selection(graph);
            bool is_init_ratio = false;
            ui ratio_count = 0;
            ui total_iter_num = 0;
            while (ratio_selection.ratioSelection(graph.getVerticesCount(),
                                                  ratio,
                                                  is_init_ratio,
                                                  is_vw,
                                                  is_dc,
                                                  ratio_o,
                                                  ratio_p,
                                                  graph.subgraph_density,
                                                  epsilon, is_res, res_width, is_map)) {
                ratio_count++;
                Graph vw_graph(false);
                if (is_dc)
                    graph.createVertexWeightedGraph(vw_graph, (ratio.first + ratio.second) / 2);
                else
                    graph.createVertexWeightedGraph(vw_graph, ratio.first / ratio.second);
                bool flag = true;
                auto verticess = new std::vector<VertexID>[1];
                std::vector<std::vector<VertexID>> vertices(2);
                bool lp_type = alloc_type == "fista"? 1: 0;
                LinearProgramming lp(false, lp_type, 0, 0, order_type);
                vw_graph.subgraph_density_lower_bound = 0;
                vw_graph.coreReduce(vw_graph, 0, true);
                lp.Init(vw_graph);
                ui T = 1;
                while(flag) {
                    red.UndirectedkCoreReduction(vw_graph, lp, true);
                    if(ext_type == "mwu") alloc.UndirectedMWUAllocation(vw_graph, lp, T, true);
                    if(ext_type == "fw") alloc.UndirectedlpAllocation(vw_graph, lp, T, true);
                    if(ext_type == "fista") alloc.UndirectedFistaAllocation(vw_graph, lp, T, false);
                    ext.UndirectedlpExactExtraction(vw_graph, lp, verticess);
                    flag = ver.UndirectedLpAppVerification(vw_graph, lp, verticess, epsilon);
                    T <<= 1;
                }
                total_iter_num += T;
                ext.directedVWApproExtraction(graph, vw_graph, vertices, verticess, ratio_o, ratio_p);
                printf("ratio_count: %d, density: %f, S/T: %lu / %lu\n", ratio_count, graph.subgraph_density, graph.vertices[0].size(), graph.vertices[1].size());
            }
            printf("%.10lf\n",graph.subgraph_density);
            rep.set_final_density(graph.subgraph_density);
            rep.set_total_iteration(total_iter_num);
        }
    } else {
//        Reduction red;
//        Allocation alloc;
//        Extraction ext;
//        Verification ver;
        if (!is_directed) {
            double l, r;
            ui T = 1;
            FlowNetwork flow;
            if(red_type == "k-core" || red_type == "stable") graph.coreReduce(graph);
            bool lp_type = alloc_type == "fista"? 1: 0;
            LinearProgramming lp(is_directed, lp_type, 0, 0, order_type);
            if(alloc_type == "fw" || alloc_type == "fista" || alloc_type == "mwu"){
                if(order_type == 0){

                }
                if(order_type == 1) graph.coreOrder(graph);
                if(order_type == 2) graph.CoreOrder(graph);
                lp.Init(graph);
            }
            lp.Init(graph);
            l = graph.subgraph_density;
            r = graph.subgraph_density_upper_bound;
            auto vertices = new std::vector<VertexID>[1];
            bool flag = true;
            while (flag) {
                std::cout<<"T = "<<T<<std::endl;
                T <<= 1;
                if (red_type == "k-core")
                    red.UndirectedkCoreReduction(graph, lp);
                if (alloc_type == "flow-exact")
                    alloc.UndirectedflowExactAllocation(graph, flow, l, r);
                if (alloc_type == "fw")
                    alloc.UndirectedlpAllocation(graph, lp, T, is_seq);
                if (alloc_type == "fista")
                    alloc.UndirectedFistaAllocation(graph, lp, T, is_seq);
                if (alloc_type == "mwu")
                    alloc.UndirectedMWUAllocation(graph, lp, T, is_seq);
                if (ext_type == "flow-exact")
                    ext.UndirectedflowExactExtraction(graph, flow, l, r, vertices);
                if (ext_type == "cp")
                    ext.UndirectedlpExactExtraction(graph, lp, vertices);
                if (ver_type == "flow-exact")
                    flag = ver.UndirectedflowExactVerification(graph, l, r);
                if (ver_type == "cp" && red_type != "stable")
                    flag = ver.UndirectedlpVerification(graph, lp, flow, vertices);
                if(ver_type == "cp" && red_type == "stable")
                    flag = ver.UndirectedlpVerification(graph, lp, flow, vertices, true);
            }
        } else {
            //The generation of Ratio set needs to be refined.
            //How to combine divide-and-conquer strategy with our current framework
            //needs to be considered.

            std::pair<double, double> ratio;
            double ratio_o = 0, ratio_p = 0;
            double reduction_ratio = 0;
            double total_edges_num = 0;
            double total_vertices_num = 0;
            double total_iter_num = 0;
            RatioSelection ratio_selection(graph);
            bool is_init_ratio = false;
            ui ratio_count = 0;
            ui red_count = 0;

            if(is_init_wcore){
                wcore_init_density(graph);
            }

            while (ratio_selection.ratioSelection(graph.getVerticesCount(),
                                                  ratio,
                                                  is_init_ratio,
                                                  is_vw,
                                                  is_dc,
                                                  ratio_o,
                                                  ratio_p,
                                                  graph.subgraph_density,
                                                  0, is_res, res_width, is_map)) {
                ratio_count++;
                bool flag = true;
                bool is_init_red = false;
                bool is_init_lp = false;
                bool is_reduced = false;
                bool is_stable_set = false;
                double rho, rho_c;
                double l = learning_rate * graph.subgraph_density;
                double r = graph.subgraph_density_upper_bound;
                FlowNetwork flow;
                int lp_type= alloc_type == "fista"? 1: 0;
                LinearProgramming lp(is_directed, lp_type, 0, 0, order_type);
                std::vector<std::pair<VertexID, VertexID>> edges;
                WCore w_core;
                std::vector<std::vector<VertexID>> vertices(2);
                std::pair<ui, ui> best_pos(0, 0);
                static Graph subgraph(is_directed, 0);
                subgraph.subgraph_density = graph.subgraph_density;
                rep.add_density(total_iter_num, graph.subgraph_density);
                if (red_type != "appro-xy-core" && red_type != "exact-xy-core"){
                    subgraph = Graph(is_directed,0);
                }
                ui iter_count = 0;

                if (is_print_c)
                    printf("%f\n", (ratio.first + ratio.second) / 2);
                if (is_wcore_shrink) {
                    static double last_density = 0;
                    if (graph.subgraph_density > last_density * 1.25) {
                        last_density = graph.subgraph_density;
                        
                        auto begin_wcore = std::chrono::steady_clock::now();
                        if(is_log)
                            printf("Enter Wcore, edges #: %d, density %.10lf\n", graph.getEdgesCount(), graph.subgraph_density);

                        Graph _wcore = Graph(is_directed, 0);
                        w_core.getWCore(graph, _wcore, ceil(graph.subgraph_density*graph.subgraph_density/4.0));
                        graph.deg_[0] = _wcore.deg_[0];
                        graph.deg_[1] = _wcore.deg_[1];
                        graph.adj_[0] = _wcore.adj_[0];
                        graph.adj_[1] = _wcore.adj_[1];
                        graph.edges_count_ = _wcore.edges_count_;
                        is_reduced = is_init_red = false;
                        
                        auto end_wcore = std::chrono::steady_clock::now();
                        if(is_log)
                            printf("Exit Wcore, edges #: %d\n", graph.getEdgesCount());

                        rep.add_graph_size(graph.getVerticesCount(), graph.getEdgesCount()); // Report
                        rep.add_total_wcore_time(std::chrono::duration<double>(end_wcore - begin_wcore).count());
                    }
                }
                while (flag) {
                    iter_count++;
                    if (!is_reduced || (alloc_type != "fw" && alloc_type != "fista")) {
                        is_reduced = true;
                        auto begin_core = std::chrono::steady_clock::now();
                        bool is_copy = alloc_type == "flow-exact";
                        if (red_type == "exact-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, is_map, true, true, is_res, res_width, false);
                        } else if (red_type == "appro-xy-core") {
                            red.xyCoreReduction(graph, subgraph, ratio, l, r, is_init_red,
                                                is_dc, is_map, false, true, is_res, res_width, is_copy);
                        } else {
                            subgraph = graph;
                        }
                        if (is_stats) {
                            red_count++;
                            total_vertices_num += subgraph.getVerticesCount();
                            total_edges_num += subgraph.getEdgesCount();
                        }
                        if (subgraph.getEdgesCount() == 0) {
                            double c;
                            if (is_map) {
                                if (ratio.first < 1 && ratio.second > 1) {
                                    c = 1;
                                } else if (ratio.second <= 1) {
                                    c = (ratio.first + ratio.second) / 2;
                                } else if (ratio.first >= 1) {
                                    c = 2 / (1 / ratio.first + 1 / ratio.second);
                                }
                            } else
                                c = (ratio.first + ratio.second) / 2;
                            if (is_res) {
                                ratio_o = std::max(ratio.first, c / res_width);
                                ratio_p = std::min(ratio.second, c * res_width);
                            } else {
                                ratio_o = ratio.first;
                                ratio_p = ratio.second;
                            }
                            break;
                        }
                        auto end_core = std::chrono::steady_clock::now();
                        rep.add_total_xycore_time(std::chrono::duration<double>(end_core - begin_core).count());
                    }
                    if (is_mul && (alloc_type == "fw"||alloc_type=="fista") && subgraph.subgraph_density < graph.subgraph_density) {
                        auto begin_red = std::chrono::steady_clock::now();
                        is_init_red = false;
                        red.xyCoreReduction(subgraph, subgraph, ratio, subgraph.subgraph_density, r, is_init_red, is_dc, is_map, false, false, is_res, res_width, true);
                        subgraph.subgraph_density = graph.subgraph_density;
                        if(is_stats)
                        {
                            printf("subgraph edges: %d/%d = %.4lf\n", subgraph.getEdgesCount(), graph.getEdgesCount(),
                                        subgraph.getEdgesCount() * 1.0 / graph.getEdgesCount());
                        }
                        auto end_red = std::chrono::steady_clock::now();
                        rep.add_total_xycore_time(std::chrono::duration<double>(end_red - begin_red).count());
                    }
                    if (is_stable_set && (alloc_type == "fw"||alloc_type=="fista") && e_stable) {
                        red.stableSetReduction(subgraph, lp, edges, is_stable_set, true);
                    }

                    auto begin_alloc = std::chrono::steady_clock::now();
                    if (alloc_type == "fw")
                        alloc.directedCPAllocation(subgraph, lp, iter_num, is_init_lp, ratio, !is_seq, is_exp, is_map);
                    if (alloc_type == "fista")
                        alloc.directedFistaAllocation(subgraph, lp, iter_num, is_init_lp, ratio, !is_seq, is_exp, is_map, is_log);
                    if (alloc_type == "flow-exact")
                        alloc.flowExactAllocation(subgraph, flow, ratio, l, r, is_dc, is_map);
                    auto end_alloc = std::chrono::steady_clock::now();
                    rep.add_total_iter_time(std::chrono::duration<double>(end_alloc - begin_alloc).count());

                    if (ext_type == "cp")
                        ext.directedCPExtraction(subgraph, lp, best_pos, vertices, ratio, ratio_o, ratio_p, rho, rho_c, is_map);
                    if (ext_type == "flow-exact")
                        ext.flowExactExtraction(graph, subgraph, ratio, flow, l, r, ratio_o, ratio_p, is_map);

                    auto begin_verify = std::chrono::steady_clock::now();
                    if (ver_type == "cp")
                        flag = ver.directedCPVerification(graph, subgraph, lp, best_pos, vertices, ratio, rho, rho_c, ratio_o, ratio_p, is_stable_set, edges, 0, false, is_map);
                    if (ver_type == "flow-exact")
                        flag = ver.flowExactVerification(graph, l, r);
                    auto end_verify = std::chrono::steady_clock::now();
                    rep.add_total_verify_time(std::chrono::duration<double>(end_verify - begin_verify).count());
                }

                if (alloc_type == "fw" || alloc_type == "fista")
                    total_iter_num += lp.cur_iter_num;
                else
                    total_iter_num += iter_count;
                rep.add_density(total_iter_num, graph.subgraph_density);
            }
            printf("ratio count: %d, density: %.4lf, S/T: %lu/%lu\n", ratio_count, graph.subgraph_density,
                   graph.vertices[0].size(), graph.vertices[1].size());
            if (is_stats)
            {
                printf ("avg vertices #: %f\navg edges #: %f\navg iterations #: %f\n", total_vertices_num / red_count, total_edges_num / red_count, total_iter_num / ratio_count);
                printf ("total subgraph percentage:\nvertex #: %.4lf edge #: %.4lf\n",
                            total_vertices_num * 1.0 / ratio_count / graph.getVerticesCount(), 
                            total_edges_num * 1.0 / ratio_count / graph.getEdgesCount());
            }
            rep.set_final_density(graph.subgraph_density);
            rep.set_total_iteration(total_iter_num);
        }

    }
    auto end = std::chrono::steady_clock::now();
    rep.set_total_run_time(std::chrono::duration<double>(end - begin).count());
    rep.print();
    return 0;
};