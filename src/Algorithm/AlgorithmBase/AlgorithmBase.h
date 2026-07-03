#ifndef __ALGORITHMBASE_H
#define __ALGORITHMBASE_H

#include "../../Network/Graph/Graph.h"
#include "../../Network/PathMethod/PathMethodBase/PathMethod.h"
#include "../../config.h"
using namespace std;


class AlgorithmBase {
protected:
    string algorithm_name;
    string experiment_label;  // 目前實驗的標籤 (例如 "request_cnt=80 Round=0")
    Graph graph;
    map<string, double> res;
    vector<double> cdf;
    vector<SDpair> requests;
    int time_limit, memory_total, request_cnt;
    map<SDpair, vector<Path>> paths;
    void update_res();
    double A, B, n, T, tao;
    double bar(double F);
    double Fswap(double Fa, double Fb);
    double t2F(double t);
    double F2t(double F);
    double pass_tao(double F);
    const vector<Path>& get_paths(int src, int dst);

    // === Routing trace (per-request CSV，供跨演算法對照) ===
    // 每個 request 寫一行到 ../data/log/routing_trace.csv。
    // 失敗的 request 傳空的 sv（chosen 欄位以 NA 填充）。
    // outcome: accepted / fail_fid / fail_mem / no_shape
    void log_routing_trace(int req_id, int src, int dst, const string& outcome,
                           const Shape_vector& sv, const vector<int>& purify_rounds,
                           double fidelity, double prob);
    // swap 樹深度（與 Shape::recursion_check 相同的 latest-切分規則）
    static int shape_tree_depth(const Shape_vector& sv, int left, int right);
public:
    AlgorithmBase(const Graph& graph, const vector<SDpair>& requests, const map<SDpair, vector<Path>>& paths);
    const map<string, double>& get_res() const;
    double get_res(string str);
    const vector<double>& get_cdf() const;
    string get_name();
    void set_experiment_label(const string& label) { experiment_label = label; }
    virtual ~AlgorithmBase();
    virtual void run() = 0;
};

#endif