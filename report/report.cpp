#include "report.h"
Report::Report(){
    graph_size.clear();
    density.clear();
    total_run_time = 0;
    final_density = 0;
    total_iteration = 0;
    total_iter_time = 0;
    total_wcore_time = 0;
    total_xycore_time = 0;
    total_verify_time = 0;
    total_edges_num = 0;
}
void Report::add_total_iter_time(double t){
    total_iter_time += t;
}
void Report::add_total_xycore_time(double t){
    total_xycore_time += t;
}
void Report::add_total_wcore_time(double t){
    total_wcore_time += t;
}
void Report::add_total_verify_time(double t){
    total_verify_time += t;
}
void Report::add_graph_size(int n, int m){
    graph_size.push_back(std::make_pair(n,m));
}
void Report::add_density(int n, double d){
    if(!density.empty() && density.back().second == d){
        return;
    }
    density.push_back(std::make_pair(n,d));
}
void Report::set_total_run_time(double t){
    total_run_time = t;
}
void Report::set_final_density(double d){
    final_density = d;
}
void Report::set_total_iteration(int i){
    total_iteration = i;
}
void Report::add_total_edges_num(long long n){
    total_edges_num += n;
}
void Report::print(){
    std::cout << "Graph size: " << std::endl;
    for(auto i: graph_size){
        std::cout << i.first << " " << i.second << std::endl;
    }
    std::cout << "Density: " << std::endl;
    for(auto i: density){
        std::cout << i.first << " " << i.second << std::endl;
    }
    printf("Total edges num: %lld\n", total_edges_num);
    printf("Details: iter: %.2lf, xycore: %.2lf, wcore: %.2lf, verify: %.2lf, other: %.2lf\n", total_iter_time, total_xycore_time, total_wcore_time, total_verify_time, total_run_time - total_iter_time - total_xycore_time - total_wcore_time - total_verify_time);
    printf("Percentage: iter: %.2lf, xycore: %.2lf, wcore: %.2lf, verify: %.2lf, other: %.2lf\n", total_iter_time / total_run_time * 100, total_xycore_time / total_run_time * 100, total_wcore_time / total_run_time * 100, total_verify_time / total_run_time * 100, (total_run_time - total_iter_time - total_xycore_time - total_wcore_time - total_verify_time) / total_run_time * 100);
    printf("Total run time: %.4lf\n", total_run_time);
    printf("Final density: %.4lf\n", final_density);
    printf("Total iteration: %d\n", total_iteration);
}
