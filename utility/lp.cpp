#include "lp.h"

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
        result = 0;
        for (ui i = 0; i < n; i++) {
            result += 1/sqrt(ratio)* r[0][i] * r[0][i] + sqrt(ratio)* r[1][i] * r[1][i];
        }
        result *= 0.25;
        w.assign(m, std::make_pair(0, 0));
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

void LinearProgramming::RACIterate(double learning_rate, double ratio, bool is_synchronous) {
    cur_iter_num++;
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
    for(auto i:perm){
        int u = alpha[i].id_first, v = alpha[i].id_second;
        beta[i].weight_first = lr2 * w[i].first + alpha[i].weight_first;
        beta[i].weight_second = lr2 * w[i].second + alpha[i].weight_second;

        r[0][u] += 2 * sqrr * lr2 * w[i].first;
        r[1][v] += 2 / sqrr * lr2 * w[i].second;

        double eta = edges_count_ * learning_rate * 2;
        auto [w1, w2] = Proj(beta[i].weight_first - 1/(2*eta) * r[0][u],
                             beta[i].weight_second - 1/(2*eta) * r[1][v]);
        double tmp = 1/lr2 - edges_count_/learning_rate;
        w[i].first -= tmp * (w1 - alpha[i].weight_first);
        w[i].second -= tmp * (w2 - alpha[i].weight_second);
        alpha[i].weight_first = w1;
        alpha[i].weight_second = w2;
    }
    last_result = result;
    r.assign(2, std::vector<double>(nodes_count_, 0));
    for (ui i = 0; i < edges_count_; i++) {
        r[0][alpha[i].id_first] += 2 * sqrr * (lr2*w[i].first+alpha[i].weight_first);
        r[1][alpha[i].id_second] += 2 / sqrr * (lr2*w[i].second+alpha[i].weight_second);
    }
    result = 0;
    for (ui i = 0; i < nodes_count_; i++) {
        result += 1/sqrr * r[0][i] * r[0][i] + sqrr * r[1][i] * r[1][i];
    }
    result *= 0.25;
    r.assign(2, std::vector<double>(nodes_count_, 0));
    for (ui i = 0; i < edges_count_; i++) {
        r[0][alpha[i].id_first] += 2 * sqrr * alpha[i].weight_first;
        r[1][alpha[i].id_second] += 2 / sqrr * alpha[i].weight_second;
    }
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