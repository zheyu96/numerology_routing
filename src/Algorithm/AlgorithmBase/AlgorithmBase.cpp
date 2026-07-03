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
    res["actual_req_cnt"] = graph.get_actual_req_cnt();
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

int AlgorithmBase::shape_tree_depth(const Shape_vector& sv, int left, int right) {
    if(left >= right - 1) return 0;  // 單條 link
    int latest = left + 1;
    for(int i = left + 1; i < right; i++) {
        if(sv[i].second[0].second > sv[latest].second[0].second)
            latest = i;
    }
    return 1 + max(shape_tree_depth(sv, left, latest),
                   shape_tree_depth(sv, latest, right));
}

vector<AlgorithmBase::TraceRow> AlgorithmBase::trace_buffer;

void AlgorithmBase::log_routing_trace(int req_id, int src, int dst, const string& outcome,
                                      const Shape_vector& sv, const vector<int>& purify_rounds,
                                      double fidelity, double prob) {
    TraceRow row;
    row.algo = algorithm_name;
    row.exp_label = experiment_label;
    row.outcome = outcome;
    row.req_id = req_id;
    row.src = src;
    row.dst = dst;
    row.sp_hop = graph.distance(src, dst);
    row.has_shape = !sv.empty();
    row.chosen_hop = row.tree_depth = row.finish_t = -1;
    row.in_path_set = false;
    row.fidelity = fidelity;
    row.prob = prob;

    if(row.has_shape) {
        row.chosen_hop = (int)sv.size() - 1;
        // 選出的節點序列是否落在預先給的候選 path set 裡（正反向都算）
        vector<int> seq;
        for(const auto& nd : sv) seq.push_back(nd.first);
        for(const Path& p : get_paths(src, dst)) {
            if(p == seq) { row.in_path_set = true; break; }
            Path rp(p.rbegin(), p.rend());
            if(rp == seq) { row.in_path_set = true; break; }
        }
        row.tree_depth = shape_tree_depth(sv, 0, (int)sv.size() - 1);
        for(int i = 0; i < (int)seq.size(); i++) {
            if(i) row.path_str += " -> ";
            row.path_str += to_string(seq[i]);
        }
        for(int li = 0; li < row.chosen_hop; li++) {
            int r = (li < (int)purify_rounds.size()) ? purify_rounds[li] : 0;
            if(li) row.pur_str += '|';
            row.pur_str += to_string(r);
        }
        row.finish_t = 0;
        for(const auto& nd : sv)
            for(const auto& rng : nd.second)
                row.finish_t = max(row.finish_t, rng.second);
    }

    #pragma omp critical(routing_trace_buffer)
    {
        trace_buffer.push_back(row);
    }
}

void AlgorithmBase::flush_routing_trace() {
    if(trace_buffer.empty()) return;

    // 按 req_id 分組（同組內保持演算法執行順序）
    map<int, vector<const TraceRow*>> by_req;
    for(const auto& row : trace_buffer)
        by_req[row.req_id].push_back(&row);

    // --- 人讀版：同一個 request 的各演算法路徑排在一起 ---
    {
        ofstream fout("../data/log/routing_trace.txt", ios::app);
        if(fout.is_open()) {
            fout << "=== Experiment: " << trace_buffer.front().exp_label << " ===" << endl;
            for(auto& [rid, rows] : by_req) {
                const TraceRow* first = rows.front();
                fout << "req#" << rid << "  SD=(" << first->src << "," << first->dst << ")"
                     << "  sp_hop=" << first->sp_hop << endl;
                for(const TraceRow* r : rows) {
                    fout << "    " << r->algo;
                    for(int pad = (int)r->algo.size(); pad < 12; pad++) fout << ' ';
                    fout << ": " << r->outcome;
                    if(r->has_shape) {
                        fout << "  | path: " << r->path_str
                             << "  pur[" << r->pur_str << "]"
                             << "  hop=" << r->chosen_hop
                             << " detour=" << (r->chosen_hop - r->sp_hop)
                             << " in_path_set=" << (r->in_path_set ? 1 : 0)
                             << " depth=" << r->tree_depth
                             << " fid=" << r->fidelity
                             << " prob=" << r->prob
                             << " finish_t=" << r->finish_t;
                    }
                    fout << endl;
                }
            }
            fout << "----------------------------------------" << endl;
        } else {
            cerr << "[Warning] Unable to open ../data/log/routing_trace.txt" << endl;
        }
    }

    // --- 機讀版 CSV（畫圖用），欄位同先前 ---
    {
        const string csv_path = "../data/log/routing_trace.csv";
        bool need_header;
        {
            ifstream fin(csv_path);
            need_header = !fin.good() || fin.peek() == ifstream::traits_type::eof();
        }
        ofstream fout(csv_path, ios::app);
        if(fout.is_open()) {
            if(need_header) {
                fout << "algo,exp_label,req_id,src,dst,sp_hop,chosen_hop,detour,"
                     << "in_path_set,tree_depth,purify_rounds,fidelity,prob,finish_t,outcome" << endl;
            }
            for(const auto& r : trace_buffer) {
                fout << r.algo << ',' << r.exp_label << ',' << r.req_id << ','
                     << r.src << ',' << r.dst << ',' << r.sp_hop << ',';
                if(!r.has_shape) {
                    fout << "NA,NA,NA,NA,NA,NA,NA,NA," << r.outcome << endl;
                } else {
                    fout << r.chosen_hop << ',' << (r.chosen_hop - r.sp_hop) << ','
                         << (r.in_path_set ? 1 : 0) << ',' << r.tree_depth << ','
                         << r.pur_str << ',' << r.fidelity << ',' << r.prob << ','
                         << r.finish_t << ',' << r.outcome << endl;
                }
            }
        }
    }

    trace_buffer.clear();
}