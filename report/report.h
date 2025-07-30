#include<iostream>
#include<vector>
class Report{
    private:
        std::vector<std::pair<int, int> > graph_size;
        std::vector<std::pair<int,double> > density; // iter, density
        double total_run_time;
        double final_density;
        int total_iteration;
        double total_iter_time;
        double total_xycore_time;
        double total_wcore_time;
        double total_verify_time;
        long long total_edges_num;
    public:
        Report();
        void add_graph_size(int n, int m);
        void add_density(int n, double d);
        void print(); 
        void set_total_run_time(double t);
        void set_final_density(double d);
        void set_total_iteration(int i);
        void add_total_iter_time(double t);
        void add_total_xycore_time(double t);
        void add_total_wcore_time(double t);
        void add_total_verify_time(double t);
        void add_total_edges_num(long long n);
};