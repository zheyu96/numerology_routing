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

    // === Routing trace (per-request，供跨演算法對照) ===
    // 各演算法把每個 request 的結果 push 進共用 buffer；同一輪所有演算法
    // 跑完後由 main 呼叫 flush_routing_trace()，按 request 分組寫出
    // routing_trace.txt（人讀，同 request 的各演算法路徑排在一起）與
    // routing_trace.csv（畫圖用）。
    // outcome: accepted / fail_fid / fail_mem / no_shape
    struct TraceRow {
        string algo, exp_label, outcome;
        int req_id, src, dst, sp_hop;
        bool has_shape;
        int chosen_hop, tree_depth, finish_t;
        bool in_path_set;
        string path_str;   // 節點序列，如 "3 -> 17 -> 45"
        string pur_str;    // 每段 purify 輪數，如 "1|0"
        double fidelity, prob;
    };
    static vector<TraceRow> trace_buffer;
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
    // 該輪所有演算法跑完後呼叫：分組寫檔並清空 buffer
    static void flush_routing_trace();
    virtual ~AlgorithmBase();
    virtual void run() = 0;
};

#endif