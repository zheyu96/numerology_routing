#include "AlgorithmBase.h"

AlgorithmBase::AlgorithmBase(const Graph& graph, const vector<SDpair>& requests, const map<SDpair, vector<Path>>& paths):
    graph(graph), requests(requests), paths(paths) {
    time_limit = this->graph.get_time_limit();
    A = this->graph.get_A();
    B = this->graph.get_B();
    n = this->graph.get_n();
    T = this->graph.get_T();
    tao = this->graph.get_tao();

    vector<bool> passed_node(this->graph.get_num_nodes(), false);
    memory_total = 0;
    // for(int i = 0; i < (int)requests.size(); i++) {
    //     int src = requests[i].first, dst = requests[i].second;
    //     vector<Path> path = graph.get_path(src, dst);
    //     for(auto node : path) {
    //         if(!passed_node[node]) {
    //             passed_node[node] = true;
    //             memory_total += graph.get_node_memory(node);
    //         }
    //     }
    // }

    memory_total *= this->graph.get_time_limit();
    request_cnt = this->requests.size();
    cdf.resize(this->graph.get_boundary().size(), request_cnt);
}
AlgorithmBase::~AlgorithmBase() {
    if(DEBUG) cerr << "delete AlgorithmBase" << endl;
}

string AlgorithmBase::get_name() {
    return algorithm_name;
}

const map<string, double>& AlgorithmBase::get_res() const {
    return res;
}

double AlgorithmBase::get_res(string str) {
    return res[str];
}
const vector<double>& AlgorithmBase::get_cdf() const {
    return cdf;
}

double AlgorithmBase::bar(double F) {
    return (1.0 - F);
}
double AlgorithmBase::Fswap(double Fa, double Fb) {
    if(Fa <= A + EPS || Fb <= A + EPS) return 0;
    return Fa * Fb + (1.0 / 3.0) * bar(Fa) * bar(Fb);
}
double AlgorithmBase::t2F(double t) {
    if(t >= INF / 1e2) return 0;
    return A + B * exp(-pow(t / T, n));
}
double AlgorithmBase::F2t(double F) {
    if(F <= A + EPS) return INF;
    return T * pow(-log((F - A) / B), 1.0 / n);
}
double AlgorithmBase::pass_tao(double F) {
    return t2F(F2t(F) + tao);
}

void AlgorithmBase::update_res() {
    res["fidelity_gain"] = graph.get_fidelity_gain();
    res["succ_request_cnt"] = graph.get_succ_request_cnt();
    res["utilization"] = (double)graph.get_usage() / (double)memory_total;
    res["pure_fidelity"] = graph.get_pure_fidelity();
    cdf.clear();
    const vector<double>& boundary = graph.get_boundary();
    const vector<double>& cnt = graph.get_cnt();
    cdf.resize(boundary.size(), 0);
    double unfinish = request_cnt;
    for(int i = 0; i < (int)boundary.size(); i++) {
        unfinish -= cnt[i];
    }
    cdf[0] = cnt[0] + unfinish;
    for(int i = 1; i < (int)boundary.size(); i++) {
        cdf[i] = cdf[i - 1] + cnt[i];
    }
}

const vector<Path>& AlgorithmBase::get_paths(int src, int dst) {
    return paths[{src, dst}];
}