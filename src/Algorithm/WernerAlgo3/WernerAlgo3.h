#ifndef __WERNER_ALGO3_H
#define __WERNER_ALGO3_H

#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"

#include <map>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <limits>
#include <cassert>

using namespace std;

/**
 * WernerAlgo3 = ZFA2 的 Upper Bound 版本（ZFA_UB）
 * - LP relaxation 階段和 ZFA2 完全相同（含 purification DP）
 * - 不做 greedy rounding，直接從 fractional solution 計算 LP 上界
 * - 關係：WernerAlgo2 (ZFA2) → WernerAlgo3 (ZFA_UB)
 *         如同 MyAlgo1 → MyAlgo2
 */
class WernerAlgo3 : public AlgorithmBase {
public:
    #define double long double
    WernerAlgo3(const Graph& graph,
               const vector<pair<int,int>>& requests,
               const map<SDpair, vector<Path>>& paths);

    void run();

private:
    // ===== Werner DP Label =====
    enum class Op : unsigned char { LEAF = 0, CONT = 1, MERGE = 2 };
    class ZLabel {
        public:
        double B= 1000.0;
        double Z = 1000.0;
        double P=1.0;
        int purify_type=-1;
        int a = -1, b = -1, t = -1, k = -1;
        int left_id=-1,right_id=-1,parent_id=-1;
        Op op = Op::LEAF;
        vector<int> ent_time;
        ZLabel(){}
        ZLabel(double _B, double _Z, double _P, Op _op,int _purify_type, int _a, int _b, int _t, int _k, int pid = -1, int lid = -1, int rid = -1)
        : B(_B), Z(_Z), P(_P), op(_op),purify_type(_purify_type), a(_a), b(_b), t(_t), k(_k), parent_id(pid), left_id(lid), right_id(rid) {}
    };

    struct DPParam{
        double eps_bucket,Zhat,Zmin,eta,T,deltaP;
        int tau_max;
    }dpp;
    double epsilon = 0.5;
    double obj = 0.0;
    vector<double> alpha;
    vector<vector<double>> beta;
    vector<map<Shape_vector,double>> x;
    vector<vector<vector<vector<ZLabel>>>> DP_table;

    void variable_initialize();
    Shape_vector separation_oracle();
    void run_dp_in_t(const Path& path, const DPParam& dpp,int t);
    void pareto_prune_byZ(vector<ZLabel>& cand);
    void bucket_by_ZP(vector<ZLabel>& cand);
    Shape_vector backtrack_shape(ZLabel leaf, const vector<int>& path, vector<int>& out_purify_rounds);
    int split_dis(int s,int d,WernerAlgo3::ZLabel& L);
    pair<double,WernerAlgo3::ZLabel> eval_best_J(int s, int d, int t, double alp);
    int purify_time=5;
    double Purify_in_vt[6][7]={
        {1,1},
        {1,2,2},
        {1,2,3,2},
        {1,2,3,3,2},
        {1,2,3,3,3,2},
        {1,2,3,3,3,3,2},
    };

    ZLabel gen_leaf_label(int s,int e,int st,int tlen,int path_a,int path_b);
    map<Shape_vector, vector<int>> shape_purify_map;

    // --- Oracle cache for incremental separation_oracle ---
    struct OracleCache {
        double best_score = 1e18;
        Shape_vector shape;
        vector<int> purify_rounds;
        bool valid = false;
    };
    vector<vector<OracleCache>> oracle_cache;
    set<int> dirty_nodes;
    set<int> dirty_alpha_idxs;
};

#endif // __WERNER_ALGO3_H
