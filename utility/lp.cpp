#include "lp.h"
#include <thread>
#include <chrono>

LinearProgramming::LinearProgramming(bool is_directed, ui type, ui n, ui m, ui sort)
        : // 0 stands for FW 1 stands for FISTA
        is_directed_(is_directed),
        nodes_count_(n),
        edges_count_(m),
        type_(type),
        sort_type(sort),
        cur_iter_num(0){
    if (is_directed_) {
        r.resize(2);
        alpha.resize(static_cast<unsigned long>(m)), 0;
        for (int i = 0; i < 2; i++) {
            r[i].resize(static_cast<unsigned long>(n), 0);
        }
    } else {
        r.resize(1);
        r[0].resize(static_cast<unsigned long>(n));
        alpha.resize(static_cast<unsigned long>(m));

        if (type_ == 1) {
            beta.resize(static_cast<unsigned long>(m));
        }
    }
}

LinearProgramming::~LinearProgramming() {

}

void LinearProgramming::InitRAC(Graph &graph, double ratio) {
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    nodes_count_ = n;
    edges_count_ = m;
    cur_iter_num = 0;
    weight = graph.weight_;
    if (is_directed_) {
        perm.resize(m);
        for(ui i = 0; i < m; i++)
            perm[i] = i;

        ui cnt = 0;
        r.assign(2, std::vector<double>(n, 0));
        sx.assign(2, std::vector<double>(n, 0));
        sy.assign(2, std::vector<double>(n, 0));
        sw.assign(2, std::vector<double>(n, 0));
        w.assign(m, std::make_pair(0, 0));
        y.resize(m);
        z.resize(m);
        for (VertexID u = 0; u < n; u++) {
            for (auto &v: graph.getOutNeighbors(u)) {
                z[cnt].weight_first = 0.5;
                z[cnt].weight_second = 0.5;
                z[cnt].id_first = u;
                z[cnt].id_second = v;
                y[cnt] = z[cnt];
                sy[0][u] += y[cnt].weight_first;
                sy[1][v] += y[cnt].weight_second;
                sx[0][u] += z[cnt].weight_first;
                sx[1][v] += z[cnt].weight_second;
                cnt++;
            }
        }

        last_result = result = 0;
        for (ui i = 0; i < n; i++) {
            result += 2 * sqrt(ratio) * sx[0][i] * sx[0][i] + 2 / sqrt(ratio) * sx[1][i] * sx[1][i];
        }
        // sort(graph);
    } else {
        throw std::runtime_error("RAC is not implemented for undirected graphs");
    }
}

void LinearProgramming::Init(Graph &graph, double ratio) {
    ui n = graph.getVerticesCount();
    ui m = graph.getEdgesCount();
    nodes_count_ = n;
    edges_count_ = m;
    cur_iter_num = 0;
    weight = graph.weight_;
    if (is_directed_) {
        last_result = 0;
        ui cnt = 0;
        for (ui i = 0; i < 2; i++)
            r[i].resize(n, 0);
        alpha.resize(m);
        for (VertexID u = 0; u < n; u++) {
            for (auto &v: graph.getOutNeighbors(u)) {
                alpha[cnt].weight_first = 0.5;
                alpha[cnt].weight_second = 0.5;
                alpha[cnt].id_first = u;
                alpha[cnt].id_second = v;
                cnt++;
            }
        }
        for (ui u = 0; u < m; u++) {
            r[0][alpha[u].id_first] += 2 * sqrt(ratio) * alpha[u].weight_first;
            r[1][alpha[u].id_second] += 2 / sqrt(ratio) * alpha[u].weight_second;
        }
        if (type_ == 1) {
            beta.resize(m);
            for (ui i = 0; i < m; i++) {
                beta[i] = alpha[i];
            }
        }
        sort(graph);
    } else {
        ui cnt = 0;
        r[0].resize(n);
        alpha.resize(m);
        for (ui u = 0; u < n; u++) {
            r[0][u] = 0;
            for (auto &v: graph.getNeighbors(u)) {
                if (v > u) continue;
                alpha[cnt].weight_first = 0.5;
                alpha[cnt].weight_second = 0.5;
                alpha[cnt].id_first = u;
                alpha[cnt].id_second = v;
                cnt++;
            }
        }
        for (ui u = 0; u < m; u++) {
            r[0][alpha[u].id_first] += alpha[u].weight_first / weight[alpha[u].id_first];

            r[0][alpha[u].id_second] += alpha[u].weight_second / weight[alpha[u].id_second];
        }
        if (type_ == 1) {
            beta.resize(m);
            for (ui i = 0; i < m; i++) {
                beta[i] = alpha[i];
            }
        }
    }
}


void LinearProgramming::Iterate(double learning_rate, double ratio, bool is_synchronous) {
    cur_iter_num++;
    if (is_directed_) {
        if (is_synchronous) {
            for (ui i = 0; i < nodes_count_; i++) {
                r[0][i] *= (1 - learning_rate);
                r[1][i] *= (1 - learning_rate);
            }
            for (ui i = 0; i < edges_count_; i++) {
                alpha[i].weight_first *= (1 - learning_rate);
                alpha[i].weight_second *= (1 - learning_rate);
                if (r[0][alpha[i].id_first] < r[1][alpha[i].id_second]) {
                    alpha[i].weight_first += learning_rate;
                    r[0][alpha[i].id_first] += 2 * sqrt(ratio) * learning_rate;
                } else {
                    alpha[i].weight_second += learning_rate;
                    r[1][alpha[i].id_second] += 2 / sqrt(ratio) * learning_rate;
                }
            }
        } else {
            std::vector<Alpha> alpha_hat(edges_count_);
            for (ui i = 0; i < edges_count_; i++) {
                alpha[i].weight_first *= (1 - learning_rate);
                alpha[i].weight_second *= (1 - learning_rate);
                if (r[0][alpha[i].id_first] < r[1][alpha[i].id_second]) {
                    alpha[i].weight_first += learning_rate;
//                    r[0][alpha[i].id_first] += 2 * sqrt(ratio) * learning_rate;
                } else {
                    alpha[i].weight_second += learning_rate;
//                    r[1][alpha[i].id_second] += 2 / sqrt(ratio) * learning_rate;
                }
            }
            r[0].assign(nodes_count_, 0);
            r[1].assign(nodes_count_, 0);

//            for (ui u = 0; u < nodes_count_; u++) {
//                r[0][u] = 0;
//                r[1][u] = 0;
//            }
//            for (ui u = 0; u < edges_count_; u++) {
//                r[0][alpha[u].id_first] = 0;
//                r[1][alpha[u].id_second] = 0;
//            }
            for (ui u = 0; u < edges_count_; u++) {
                r[0][alpha[u].id_first] += 2 * sqrt(ratio) * alpha[u].weight_first;
                r[1][alpha[u].id_second] += 2 / sqrt(ratio) * alpha[u].weight_second;
            }
        }
    } else {
        std::vector<Alpha> alpha_hat;
        alpha_hat.resize(edges_count_);
        if(is_synchronous){
            for (ui i = 0; i < nodes_count_; i++) r[0][i] *= (1 - learning_rate);
            for (ui i = 0; i < edges_count_; i++) {
                alpha[i].weight_first *= (1 - learning_rate);
                alpha[i].weight_second *= (1 - learning_rate);
                if (r[0][alpha[i].id_first] < r[0][alpha[i].id_second])
                    alpha[i].weight_first += learning_rate, r[0][alpha[i].id_first] += learning_rate / weight[alpha[i].id_first];
                else
                    alpha[i].weight_second += learning_rate, r[0][alpha[i].id_second] += learning_rate / weight[alpha[i].id_second];
            }
        }
        else{
            for (ui i = 0; i < edges_count_; i++) {
                if (r[0][alpha[i].id_first] < r[0][alpha[i].id_second])
                    alpha_hat[i].weight_first = 1, alpha_hat[i].weight_second = 0;
                else
                    alpha_hat[i].weight_second = 1, alpha_hat[i].weight_first = 0;
            }
            for (ui i = 0; i < edges_count_; i++) {

                alpha[i].weight_first = (1 - learning_rate) * alpha[i].weight_first
                                        + learning_rate * alpha_hat[i].weight_first;
                alpha[i].weight_second = (1 - learning_rate) * alpha[i].weight_second
                                        + learning_rate * alpha_hat[i].weight_second;
            }
            for (ui i = 0; i < nodes_count_; i++) {
                r[0][i] = 0;
            }
            for (ui i = 0; i < edges_count_; i++) {
                r[0][alpha[i].id_first] += alpha[i].weight_first / weight[alpha[i].id_first];
                r[0][alpha[i].id_second] += alpha[i].weight_second / weight[alpha[i].id_second];
            }
        }
        
    }
}

void LinearProgramming::FistaIterate(double learning_rate, double t, double ratio, bool is_synchronous) {
    if (is_directed_) {
        auto Proj = [](double x,double y) -> std::pair<double,double>{
            if (fabs(x - y) <= 1) return std::make_pair((x - y + 1) / 2, (y - x + 1) / 2);
            if (x - y > 0) return std::make_pair(1.0, 0.0);
            return std::make_pair(0.0, 1.0);
        };
        ++cur_iter_num;
        std::vector<Alpha> alpha_new(edges_count_);
        double gamma_t = (t - 1) / (t + 2);
        if (is_synchronous){
            throw std::runtime_error("FISTA's synchronous mode is unverified.");
            // for (ui i = 0; i < edges_count_; i++) {
            //     beta[i].weight_first = beta[i].weight_first - learning_rate * r[0][beta[i].id_first];
            //     beta[i].weight_second = beta[i].weight_second - learning_rate * r[1][beta[i].id_second];

            //     beta[i].weight_first += learning_rate * ratio;
            //     beta[i].weight_second += learning_rate / ratio;

            //     auto [w1, w2] = Proj(beta[i].weight_first, beta[i].weight_second);
            //     beta[i].weight_first = w1;
            //     beta[i].weight_second = w2;
                
            //     r[0][beta[i].id_first] -= 2 * sqrt(ratio) * alpha[i].weight_first;
            //     r[1][beta[i].id_second] -= 2 / sqrt(ratio) * alpha[i].weight_second;
            //     r[0][beta[i].id_first] += 2 * sqrt(ratio) * beta[i].weight_first;
            //     r[1][beta[i].id_second] += 2 / sqrt(ratio) * beta[i].weight_second;

            //     alpha_new[i] = beta[i];
            //     beta[i].weight_first =  alpha_new[i].weight_first + (alpha_new[i].weight_first - alpha[i].weight_first) * gamma_t;
            //     beta[i].weight_second = alpha_new[i].weight_second + (alpha_new[i].weight_second - alpha[i].weight_second) * gamma_t;
            // }
            // alpha = alpha_new;
        }
        else{
            for (ui i = 0; i < edges_count_; i++) {
                beta[i].weight_first = beta[i].weight_first - learning_rate * r[0][beta[i].id_first];
                beta[i].weight_second = beta[i].weight_second - learning_rate * r[1][beta[i].id_second];

                auto [w1, w2] = Proj(beta[i].weight_first, beta[i].weight_second);
                beta[i].weight_first = w1;
                beta[i].weight_second = w2;
            }
            
            alpha_new = beta;
            for (ui i = 0; i < edges_count_; i++) {
                beta[i].weight_first = alpha_new[i].weight_first + (alpha_new[i].weight_first - alpha[i].weight_first) * gamma_t;
                beta[i].weight_second = alpha_new[i].weight_second + (alpha_new[i].weight_second - alpha[i].weight_second) * gamma_t;
            }
            alpha = alpha_new;
            r[0].assign(nodes_count_, 0);
            r[1].assign(nodes_count_, 0);
            for (ui i = 0; i < edges_count_; i++) {
                r[0][alpha[i].id_first] += beta[i].weight_first;
                r[1][alpha[i].id_second] += beta[i].weight_second;
            }
            for (ui i = 0; i < nodes_count_; i++) {
                r[0][i] = r[0][i] * 2 * sqrt(ratio);
                r[1][i] = r[1][i] * 2 / sqrt(ratio);
            }
        }
    } else {
        if(is_synchronous){
            std::vector<Alpha> alpha_new;
            alpha_new.resize(edges_count_);
            for (ui i = 0; i < edges_count_; i++) {
                beta[i].weight_first = beta[i].weight_first - 2 * learning_rate * r[0][beta[i].id_first];
                beta[i].weight_second = beta[i].weight_second - 2 * learning_rate * r[0][beta[i].id_second];
                if (abs(beta[i].weight_first - beta[i].weight_second) < 1) {
                    beta[i].weight_first = (beta[i].weight_first - beta[i].weight_second + 1) / 2;
                    beta[i].weight_second = 1 - beta[i].weight_first;
                } else if (beta[i].weight_first - beta[i].weight_second > 0) {
                    beta[i].weight_first = 1;
                    beta[i].weight_second = 0;
                } else {
                    beta[i].weight_first = 0;
                    beta[i].weight_second = 1;
                }
                r[0][beta[i].id_first] -= alpha[i].weight_first / weight[beta[i].id_first];
                r[0][beta[i].id_second] -= alpha[i].weight_second / weight[beta[i].id_second];
                r[0][beta[i].id_first] += beta[i].weight_first / weight[beta[i].id_first];
                r[0][beta[i].id_second] += beta[i].weight_second / weight[beta[i].id_second];
                alpha_new[i] = beta[i];
                beta[i].weight_first =
                        alpha_new[i].weight_first + (alpha_new[i].weight_first - alpha[i].weight_first) * (t - 1) / (t + 2);
                beta[i].weight_second = alpha_new[i].weight_second +
                                        (alpha_new[i].weight_second - alpha[i].weight_second) * (t - 1) / (t + 2);
            }  
            alpha = alpha_new;
        }
        else{
            std::vector<Alpha> alpha_new;
            for (ui i = 0; i < edges_count_; i++) {
                beta[i].weight_first = beta[i].weight_first - 2 * learning_rate * r[0][beta[i].id_first];
                beta[i].weight_second = beta[i].weight_second - 2 * learning_rate * r[0][beta[i].id_second];
                if (abs(beta[i].weight_first - beta[i].weight_second) < 1) {
                    beta[i].weight_first = (beta[i].weight_first - beta[i].weight_second + 1) / 2;
                    beta[i].weight_second = 1 - beta[i].weight_first;
                } else if (beta[i].weight_first - beta[i].weight_second > 0) {
                    beta[i].weight_first = 1;
                    beta[i].weight_second = 0;
                } else {
                    beta[i].weight_first = 0;
                    beta[i].weight_second = 1;
                }
            }
            alpha_new = beta;
            for (ui i = 0; i < edges_count_; i++) {
                beta[i].weight_first =
                        alpha_new[i].weight_first + (alpha_new[i].weight_first - alpha[i].weight_first) * (t - 1) / (t + 2);
                beta[i].weight_second = alpha_new[i].weight_second +
                                        (alpha_new[i].weight_second - alpha[i].weight_second) * (t - 1) / (t + 2);
            }
            alpha = alpha_new;
            for (ui i = 0; i < nodes_count_; i++) {
                r[0][i] = 0;
            }
            for (ui i = 0; i < edges_count_; i++) {
                r[0][alpha[i].id_first] += alpha[i].weight_first / weight[alpha[i].id_first];
                r[0][alpha[i].id_second] += alpha[i].weight_second / weight[alpha[i].id_second];
            }
        }
    }
}

void LinearProgramming::RACRestart(double &learning_rate, double ratio) {
    sw.assign(2, std::vector<double>(nodes_count_, 0));
    sx.assign(2, std::vector<double>(nodes_count_, 0));
    sy.assign(2, std::vector<double>(nodes_count_, 0));
    double lr2= learning_rate * learning_rate;
    for(int i = 0; i < edges_count_; i++){
        z[i].weight_first = lr2*w[i].first + y[i].weight_first;
        z[i].weight_second = lr2*w[i].second + y[i].weight_second;
        y[i] = z[i];
        w[i] = std::make_pair(0,0);
        int u = y[i].id_first, v = y[i].id_second;
        sy[0][u] += y[i].weight_first;
        sy[1][v] += y[i].weight_second;
        sx[0][u] += y[i].weight_first;
        sx[1][v] += y[i].weight_second;
        sw[0][u] += w[i].first;
        sw[1][v] += w[i].second;
    }
}

void LinearProgramming::RACIterate(double learning_rate, double ratio, bool is_synchronous) {
    if (!is_directed_){
        throw std::runtime_error("RAC is not implemented for undirected graphs.");
    }
    if (is_synchronous) {
        throw std::runtime_error("RAC's synchronous mode is unverified.");
    }
    ++cur_iter_num;
    auto Proj = [](double x, double y) -> std::pair<double, double> {
        if (fabs(x - y) <= 1) return std::make_pair((x - y + 1) / 2, (y - x + 1) / 2);
        if (x - y > 0) return std::make_pair(1.0, 0.0);
        return std::make_pair(0.0, 1.0);
    };
    double sqrr = sqrt(ratio);
    double lr2 = learning_rate * learning_rate;
    ui m = edges_count_;
    for(auto i:perm){ 
        int u = y[i].id_first, v = y[i].id_second;
        double eta = 2 * m * learning_rate * std::max(sqrr,1/sqrr);
        auto [y1, y2] = Proj(y[i].weight_first  - 1/(2*eta)*2*sqrr*sx[0][u],
                             y[i].weight_second - 1/(2*eta)*2/sqrr*sx[1][v]);
        
        double w1 = w[i].first  - (1 - m * learning_rate) / lr2 * (y1 - y[i].weight_first);
        double w2 = w[i].second - (1 - m * learning_rate) / lr2 * (y2 - y[i].weight_second);

        last_result = result;
        
        result -= 2 * sqrr * sx[0][u] * sx[0][u] + 2 / sqrr * sx[1][v] * sx[1][v];

        sw[0][u] = sw[0][u] + w1 - w[i].first;
        sw[1][v] = sw[1][v] + w2 - w[i].second;
        sy[0][u] = sy[0][u] + y1 - y[i].weight_first;
        sy[1][v] = sy[1][v] + y2 - y[i].weight_second;

        sx[0][u] = lr2 * sw[0][u] + sy[0][u];
        sx[1][v] = lr2 * sw[1][v] + sy[1][v];

        result += 2 * sqrr * sx[0][u] * sx[0][u] + 2 / sqrr * sx[1][v] * sx[1][v];

        y[i].weight_first = y1;
        y[i].weight_second = y2;
        w[i].first = w1;
        w[i].second = w2;
    }
    for (ui i = 0; i < m; ++i){
        z[i].weight_first  = lr2*w[i].first  + y[i].weight_first;
        z[i].weight_second = lr2*w[i].second + y[i].weight_second;
    }
    for (ui i = 0; i < nodes_count_; ++i){
        r[0][i] = 2 * sqrr * sx[0][i];
        r[1][i] = 2 / sqrr * sx[1][i];
    }
    // -- diagnostic output
    // std::cout << "Iteration: " << cur_iter_num << ", Result: " << result << ", Last Result: " << last_result << std::endl;
    // for (ui i = 0; i < nodes_count_; ++i){
    //     std::cout << "Node " << i << ": r[0] = " << r[0][i] << ", r[1] = " << r[1][i] << std::endl;
    // }
    // for (ui i = 0; i < m; ++i){
    //     printf("Edge %d: (%d, %d) -> (%.5lf, %.5lf)\n", i, y[i].id_first, y[i].id_second, y[i].weight_first, y[i].weight_second);
    //     // printf("                  -> (%.5lf, %.5lf)\n", i, y[i].id_first, y[i].id_second, z[i].weight_first, z[i].weight_second);
    // }
    // let the program sleep for 0.5s
    // std::this_thread::sleep_for(std::chrono::milliseconds(500));
}

void LinearProgramming::MWUIterate(ui t, bool is_synchronous){
    double learning_rate = 1.0 / (t + 1);
    if(is_synchronous){
        for (ui i = 0; i < nodes_count_; i++) r[0][i] *= (1 - learning_rate);
        for (ui i = 0; i < edges_count_; i++) {
            alpha[i].weight_first *= (1 - learning_rate);
            alpha[i].weight_second *= (1 - learning_rate);
            if (r[0][alpha[i].id_first] < r[0][alpha[i].id_second])
                alpha[i].weight_first += learning_rate, r[0][alpha[i].id_first] += learning_rate / weight[alpha[i].id_first];
            else
                alpha[i].weight_second += learning_rate, r[0][alpha[i].id_second] += learning_rate / weight[alpha[i].id_second];
        }
    }
    else{
        std::vector<Alpha> alpha_hat;
        alpha_hat.resize(edges_count_);
        for (ui i = 0; i < edges_count_; i++) {
            if (r[0][alpha[i].id_first] < r[0][alpha[i].id_second])
                alpha_hat[i].weight_first = 1, alpha_hat[i].weight_second = 0;
            else
                alpha_hat[i].weight_second = 1, alpha_hat[i].weight_first = 0;
        }
        for (ui i = 0; i < edges_count_; i++) {

            alpha[i].weight_first = (1 - learning_rate) * alpha[i].weight_first
                                    + learning_rate * alpha_hat[i].weight_first;
            alpha[i].weight_second = (1 - learning_rate) * alpha[i].weight_second
                                        + learning_rate * alpha_hat[i].weight_second;
        }
        for (ui i = 0; i < nodes_count_; i++) {
            r[0][i] = 0;
        }
        for (ui i = 0; i < edges_count_; i++) {
            r[0][alpha[i].id_first] += alpha[i].weight_first / weight[alpha[i].id_first];
            r[0][alpha[i].id_second] += alpha[i].weight_second / weight[alpha[i].id_second];
        }
    }
}


void LinearProgramming::sort(Graph &graph) {
    if (is_directed_) {
        if (sort_type == 1) {
            std::vector<std::vector<ui>> degrees(2);
            degrees[0] = graph.getOutDegrees();
            degrees[1] = graph.getInDegrees();
            std::sort(alpha.begin(), alpha.end(), [this, &degrees](Alpha a, Alpha b) -> bool {
                return degrees[0][a.id_first] * degrees[1][a.id_second] <
                       degrees[0][b.id_first] * degrees[1][b.id_second];
            });
        } else if (sort_type == 2) {
            std::vector<std::vector<ui>> degrees(2);
            degrees[0] = graph.getOutDegrees();
            degrees[1] = graph.getInDegrees();
            std::sort(alpha.begin(), alpha.end(), [this, &degrees](Alpha a, Alpha b) -> bool {
                return degrees[0][a.id_first] + degrees[1][a.id_second] <
                       degrees[0][b.id_first] + degrees[1][b.id_second];
            });
        }
    } else {

    }
}