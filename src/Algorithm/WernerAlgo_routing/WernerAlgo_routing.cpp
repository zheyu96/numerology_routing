#include "../../Network/Purification/Purification.h"
#include "WernerAlgo_routing.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <climits>

using namespace std;

WernerAlgo_routing::WernerAlgo_routing(const Graph& graph,const vector<pair<int,int>>& requests,const map<SDpair, vector<Path>>& paths,double epsilon,double bucket_eps): AlgorithmBase(graph, requests, paths)
{
    algorithm_name = "ZFA_routing";
    this->epsilon = epsilon;
    this->bucket_eps = bucket_eps;
}

void WernerAlgo_routing::initialize_bucket_minima() {
    dpp.Zmin = INF;
    dpp.Pmin = INF;
    const int max_rounds = min(purify_time, graph.get_time_limit() - 3);

    for (int u = 0; u < graph.get_num_nodes(); ++u) {
        for (int v : graph.adj_list[u]) {
            if (u >= v) continue;
            const double w_ini = graph.get_link_werner(u, v);
            const double entangle_prob = graph.get_entangle_succ_prob(u, v);
            for (int rounds = 0; rounds <= max_rounds; ++rounds) {
                const double w = Purification::pumping_werner(w_ini, rounds);
                const double probability = Purification::pumping_success_prob(
                    entangle_prob, w_ini, rounds);
                if (w <= 0.0 || w > 1.0 || probability <= 0.0 || probability > 1.0)
                    continue;

                const double W = sqrt(max((double)0.0, -log(w)));
                const double P = -log(probability);
                if (W > 0.0 && W <= dpp.Zhat + EPS) dpp.Zmin = min(dpp.Zmin, W);
                if (P > 0.0) dpp.Pmin = min(dpp.Pmin, P);
            }
        }
    }

    if (dpp.Zmin >= INF / 2) dpp.Zmin = max((double)graph.get_Zmin(), (double)1e-12);
    if (dpp.Pmin >= INF / 2) dpp.Pmin = (double)1e-12;
}

void WernerAlgo_routing::variable_initialize() {
    // 與 MyAlgo1 類似：初始化 dual 與目標
    int m = (int)requests.size()
          + graph.get_num_nodes() * graph.get_time_limit();

    double delta = (1 + epsilon) * (1.0 / pow((1 + epsilon) * m, 1.0 / epsilon));
    obj = m * delta;

    alpha.assign(requests.size(), delta);
    x.clear();
    x.resize(requests.size());
    int V = graph.get_num_nodes();
    int T = graph.get_time_limit();
    dpp.eps_bucket = (bucket_eps > 0.0) ? bucket_eps : graph.get_bucket_eps();
    double F_th=graph.get_fidelity_threshold();
    double w_th=(4.0*F_th-1.0)/3.0;
    dpp.Zhat = sqrt(-log(w_th))+1e-9;
    dpp.T    = time_limit-1;
    dpp.tau_max=min(time_limit-1,5);
    dpp.eta  = graph.get_tao()/graph.get_T();   // 論文 Eq.(2): eta = delta / T_mem
    initialize_bucket_minima();
    beta.assign(V, vector<double>(T, INF));

    for (int v = 0; v < V; ++v) {
        for (int t = 0; t < T; ++t) {
            int cap = graph.get_node_memory_at(v, t);
            beta[v][t] = (cap == 0) ? INF : (delta / cap);
        }
    }

    // 全點對 hop 距離（BFS from every node），給 split_dis 的平衡度 tie-break 用
    hop_dist.assign(V, vector<int>(V, INT_MAX / 2));
    for (int s = 0; s < V; ++s) {
        hop_dist[s][s] = 0;
        queue<int> que;
        que.push(s);
        while (!que.empty()) {
            int u = que.front(); que.pop();
            for (int v2 : graph.adj_list[u]) {
                if (hop_dist[s][v2] > hop_dist[s][u] + 1) {
                    hop_dist[s][v2] = hop_dist[s][u] + 1;
                    que.push(v2);
                }
            }
        }
    }

    // DP 表只在這裡配置一次，separation_oracle 內僅清空 cell 重用 buffer
    int Tcap = (int)dpp.T + 5;
    DP_table.assign(Tcap,
        vector<vector<vector<ZLabel>>>(V, vector<vector<ZLabel>>(V)));
}

Shape_vector WernerAlgo_routing::separation_oracle(){
    double most_violate=1e9;
    Shape_vector todo_shape;
    vector<int> best_purify_rounds;

    // ===== 全圖 all-pairs Werner DP（時間展開，Floyd-Warshall 風格）=====
    // DP_table[t][u][v]：節點 u,v 之間在 time t 完成的非支配候選集。
    // 注意：separation_oracle 在 run() 的 LP 迴圈中被重複呼叫，且每次呼叫前
    //   beta 都已被上一輪灌流量改過，所以這張表「每次呼叫都要用當下 beta 重算」。
    //   這裡的「共用」只指：在「同一次呼叫內」DP 只看 beta、與 request/alpha 無關，
    //   因此所有 request 共用這一張表（只算一次，再各自讀 [src][dst]）。
    // buffer 在 variable_initialize 配置一次，這裡只清空 cell（保留 capacity）。
    for(auto& layer : DP_table)
        for(auto& row : layer)
            for(auto& cell : row)
                cell.clear();
    for(int t=1;t<=dpp.T;t++)
        run_dp_in_t(dpp,t);

    // 對每個 request 直接讀 DP_table[t][src][dst]，挑最 violate 的 shape
    for(int i=0;i<(int)requests.size();i++){
        int src=requests[i].first,dst=requests[i].second;

        double local_best_J = 1e18;
        ZLabel local_best_label;
        for(int t=1;t<=dpp.T;t++){
            auto cur_val=eval_best_J(src,dst,t,alpha[i]);
            if(cur_val.first < local_best_J){
                local_best_J = cur_val.first;
                local_best_label = cur_val.second;
            }
        }

        if(local_best_J < 1e18 && local_best_J < most_violate){
            vector<int> cur_rounds;
            Shape_vector sh = backtrack_shape(local_best_label, cur_rounds);
            // all-pairs DP 可能 merge 出重複經過同一節點的 walk，這種 shape
            // 過不了 check_valid / check_resource，直接擋在 oracle 這層
            if(!sh.empty() && !has_duplicate_nodes(sh)){
                most_violate = local_best_J;
                todo_shape = sh;
                best_purify_rounds = cur_rounds;
            }
        }
    }

    if(!todo_shape.empty()){
        auto it = shape_purify_map.find(todo_shape);
        if(it == shape_purify_map.end()){
            shape_purify_map[todo_shape] = best_purify_rounds;
        } else {
            // 同一個 Shape_vector（節點+區間）可能同時來自「pur=0 leaf 經 CONT
            // 延長」與「purify leaf」：區間完全相同、key 相同。兩者都通過了
            // Zhat（預測 fidelity 皆達標），此時輪數少的 prob 高又省記憶體，
            // 對 fidelity_gain = W(fid)×pr 嚴格較優 → 保留總輪數少的版本。
            // （舊版只往 purify 方向覆寫，是單向棘輪：LP 灌在 no-purify 上的
            //   流量會在 rounding 時整批被升級成 purify 執行。）
            int new_sum = 0, old_sum = 0;
            for(int r : best_purify_rounds) new_sum += max(r, 0);
            for(int r : it->second) old_sum += max(r, 0);
            if(new_sum < old_sum)
                it->second = best_purify_rounds;
        }
    }
    return todo_shape;
}
WernerAlgo_routing::ZLabel WernerAlgo_routing::gen_leaf_label(int s,int e,int st,int tlen,int path_a,int path_b) {
    double Bleaf=0.0;
    if(st-tlen<0) return ZLabel();
    for(int i=0;i<=tlen;i++){
        double bt=beta[s][st-i]+beta[e][st-i];
        Bleaf+=bt*Purify_in_vt[tlen-1][i];
    }
    double w_ini=graph.get_link_werner(s,e);
    int rounds = tlen - 1;
    double w_cur = Purification::pumping_werner(w_ini, rounds);
    double p_cur = Purification::pumping_success_prob(graph.get_entangle_succ_prob(s,e), w_ini, rounds);
    double Zleaf=sqrt(-log(w_cur));
    double Pleaf=log(p_cur);
    if(Zleaf>dpp.Zhat) return ZLabel();
    return ZLabel(Bleaf,Zleaf,Pleaf,Op::LEAF,tlen-1,path_a,path_b,st,-1);
} 
#if 0
// Legacy capped implementation retained temporarily for comparison.
void WernerAlgo_routing::run_dp_in_t_legacy(const DPParam& dpp,int t) {
    const int n = (int)graph.get_num_nodes();

    // (t-1) 層每個 cell 的最小 Z label 位置：merge 的「可行性保底名額」。
    // 前綴改按 J 代理值 (Z²−P) 排序後，高 tau 時前綴內可能全是接近 Zhat 的
    // 低純化 label，往上 merge 全數超標；保底名額讓 purify 最兇（Z 最小）的
    // label 一定能參加 merge，oracle 才不會在長路徑上空手而回。
    vector<vector<int>> zmin_prev(n, vector<int>(n, -1));
    for(int a=0;a<n;a++)
        for(int b=0;b<n;b++){
            const auto& L = DP_table[t-1][a][b];
            for(int i=0;i<(int)L.size();i++)
                if(zmin_prev[a][b] < 0 || L[i].Z < L[zmin_prev[a][b]].Z)
                    zmin_prev[a][b] = i;
        }
    vector<int> ids1, ids2;  // merge 用的 label index buffer（迴圈間重用）

    // -------- 全圖 all-pairs：只算 a<b，[b][a] 由鏡像複製 --------
    // 網路無向、merge/leaf 的 B/Z/P 公式皆對稱，所以 [b][a] 恆為 [a][b] 的
    // 逐元素鏡像（a/b 互換、left_id/right_id 互換），維持這個不變量後
    // merge 查 [t-1][a][k]、[t-1][k][b] 不論 id 大小都直接可用，計算量減半。
    for(int a=0;a<n;a++)
        for(int b=a+1;b<n;b++){
            vector<ZLabel> cand;
            //leaf
            if(graph.adj_set[a].count(b)){
                for(int i=0;i<=purify_time;i++){
                    if(t-i-1<=0) continue;
                    ZLabel L=gen_leaf_label(a,b,t,i+1,a,b);
                    if(L.Z<=dpp.Zhat){
                        L.ent_l=t-i-1;
                        L.ent_r=t;
                        cand.push_back(L);
                    }
                }
            }
            //continue
            const auto& pre=DP_table[t-1][a][b];
            for(int p_id=0;p_id<(int)pre.size();p_id++){
                double Zp=pre[p_id].Z+dpp.eta;
                if(Zp<=dpp.Zhat){
                    double Bp=pre[p_id].B+beta[a][t]+beta[b][t];
                    double Pp=pre[p_id].P;
                    ZLabel L(Bp,Zp,Pp,Op::CONT,-1,a,b,t,-1,p_id);
                    cand.push_back(L);
                }
            }
            //merge：每個 swap 節點 k 有自己的候選配額，避免舊版「候選額度被
            //  低編號 k 先搶光、高編號節點永遠當不了 swap 點」的系統性偏差。
            //  cell 內 label 經 bucket_by_ZP 後按 J 代理值 (Z²−P) 遞增排序，
            //  取前綴即是各 k 中「eval_best_J 最可能選中」的一批；另補一個
            //  Z 最小 label 的保底名額（見 zmin_prev 註解）。
            int valid_k = 0;
            for(int k=0;k<n;k++){
                if(k==a||k==b) continue;
                if(!DP_table[t-1][a][k].empty() && !DP_table[t-1][k][b].empty())
                    valid_k++;
            }
            if(valid_k > 0){
                size_t per_k_budget = max((size_t)36, MAX_CANDIDATES_PER_CELL / (size_t)valid_k);
                int side = (int)ceil(sqrt((double)per_k_budget));
                for(int k=0;k<n;k++){
                    if(k==a||k==b) continue; //不能有自環
                    const auto& L1=DP_table[t-1][a][k];
                    const auto& L2=DP_table[t-1][k][b];
                    if(L1.size()==0||L2.size()==0) continue;
                    double swap_prob=log(graph.get_node_swap_prob(k));
                    int lim1=min((int)L1.size(), side);
                    int lim2=min((int)L2.size(), side);
                    ids1.clear(); ids2.clear();
                    for(int i=0;i<lim1;i++) ids1.push_back(i);
                    if(zmin_prev[a][k] >= lim1) ids1.push_back(zmin_prev[a][k]);
                    for(int i=0;i<lim2;i++) ids2.push_back(i);
                    if(zmin_prev[k][b] >= lim2) ids2.push_back(zmin_prev[k][b]);
                    for(int lid : ids1)
                        for(int rid : ids2){
                            const auto& left_seg=L1[lid];
                            const auto& right_seg=L2[rid];
                            double Zp=sqrt((left_seg.Z+dpp.eta)*(left_seg.Z+dpp.eta)+
                                            (right_seg.Z+dpp.eta)*(right_seg.Z+dpp.eta));
                            if(Zp<=dpp.Zhat){
                                double Pp=left_seg.P+right_seg.P+swap_prob;
                                double Bp=left_seg.B+right_seg.B+beta[a][t]+beta[b][t];
                                ZLabel L(Bp,Zp,Pp,Op::MERGE,-1,a,b,t,k,-1,lid,rid);
                                cand.push_back(L);
                            }
                        }
                }
            }
            vector<ZLabel> non_leaf;
            for (auto& L : cand) {
                if (L.op != Op::LEAF)
                    non_leaf.push_back(L);
            }
            bucket_by_ZP(non_leaf);// trimming non-leaf
            trim_cell(non_leaf);
            cand.erase(
                remove_if(cand.begin(), cand.end(),
                        [](const ZLabel& L){ return L.op != Op::LEAF; }),
                cand.end());
            cand.insert(cand.end(), non_leaf.begin(), non_leaf.end());
            trim_cell(cand);

            // 鏡像複製到 [b][a]：a/b 互換、left/right 子段索引互換，
            // 其餘欄位（含 parent_id / ent 區間 / B / Z / P）不變。
            vector<ZLabel> mirrored(cand.size());
            for(size_t i=0;i<cand.size();i++){
                ZLabel M = cand[i];
                std::swap(M.a, M.b);
                std::swap(M.left_id, M.right_id);
                mirrored[i] = M;
            }
            DP_table[t][a][b] = std::move(cand);
            DP_table[t][b][a] = std::move(mirrored);
        }
}

#endif

void WernerAlgo_routing::run_dp_in_t(const DPParam& dpp, int t) {
    const int node_count = (int)graph.get_num_nodes();

    for (int a = 0; a < node_count; ++a) {
        for (int b = a + 1; b < node_count; ++b) {
            vector<ZLabel> labels;
            map<pair<long long,long long>, ZLabel> non_base_buckets;

            auto retain_non_base = [&](const ZLabel& label) {
                const auto key = bucket_key(label);
                auto it = non_base_buckets.find(key);
                if (it == non_base_buckets.end() || label.B < it->second.B)
                    non_base_buckets[key] = label;
            };

            // E(t,[a,b]): paper base labels are kept without trimming.
            if (graph.adj_set[a].count(b)) {
                for (int rounds = 0; rounds <= purify_time; ++rounds) {
                    if (t - rounds - 1 <= 0) continue;
                    ZLabel leaf = gen_leaf_label(a, b, t, rounds + 1, a, b);
                    if (leaf.Z > dpp.Zhat) continue;
                    leaf.ent_l = t - rounds - 1;
                    leaf.ent_r = t;
                    labels.push_back(leaf);
                }
            }

            const auto& previous = DP_table[t - 1][a][b];
            for (int parent_id = 0; parent_id < (int)previous.size(); ++parent_id) {
                const double Z = previous[parent_id].Z + dpp.eta;
                if (Z > dpp.Zhat) continue;
                retain_non_base(ZLabel(
                    previous[parent_id].B + beta[a][t] + beta[b][t],
                    Z, previous[parent_id].P, Op::CONT, -1,
                    a, b, t, -1, parent_id));
            }

            // Enumerate all merge pairs and retain min-B representatives online.
            // This produces exactly the same buckets as trimming their full union.
            for (int k = 0; k < node_count; ++k) {
                if (k == a || k == b) continue;
                const auto& left = DP_table[t - 1][a][k];
                const auto& right = DP_table[t - 1][k][b];
                if (left.empty() || right.empty()) continue;
                const double swap_log_probability = log(graph.get_node_swap_prob(k));

                for (int left_id = 0; left_id < (int)left.size(); ++left_id) {
                    for (int right_id = 0; right_id < (int)right.size(); ++right_id) {
                        const double left_W = left[left_id].Z + dpp.eta;
                        const double right_W = right[right_id].Z + dpp.eta;
                        const double Z = sqrt(left_W * left_W + right_W * right_W);
                        if (Z > dpp.Zhat) continue;
                        retain_non_base(ZLabel(
                            left[left_id].B + right[right_id].B + beta[a][t] + beta[b][t],
                            Z, left[left_id].P + right[right_id].P + swap_log_probability,
                            Op::MERGE, -1, a, b, t, k, -1, left_id, right_id));
                    }
                }
            }

            for (const auto& entry : non_base_buckets)
                labels.push_back(entry.second);
            sort(labels.begin(), labels.end(), [](const ZLabel& x, const ZLabel& y) {
                return x.Z * x.Z - x.P < y.Z * y.Z - y.P;
            });

            vector<ZLabel> mirrored(labels.size());
            for (size_t i = 0; i < labels.size(); ++i) {
                mirrored[i] = labels[i];
                swap(mirrored[i].a, mirrored[i].b);
                swap(mirrored[i].left_id, mirrored[i].right_id);
            }
            DP_table[t][a][b] = std::move(labels);
            DP_table[t][b][a] = std::move(mirrored);
        }
    }
}

void WernerAlgo_routing::pareto_prune_byZ(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    sort(cand.begin(), cand.end(), [](const ZLabel& x, const ZLabel& y){
        if(x.Z!=y.Z) return x.Z < y.Z;
        return x.B<y.B;
    });
    vector<ZLabel> kept;
    double bestB = INF;
    for (auto& L : cand) {
        if (L.B + 1e-12 < bestB) {
            kept.push_back(L);
            bestB = L.B;
        }
    }
    cand.swap(kept);
}

pair<long long,long long> WernerAlgo_routing::bucket_key(const ZLabel& label) const {
    double q = 1 + dpp.eps_bucket;
    if (q <= 1.0) q = 1.0 + 1e-12;
    const double inv_log_q = 1.0 / log(q);

    long long kW = 0;
    if (label.Z > dpp.Zmin)
        kW = max(0LL, (long long)floor(log(label.Z / dpp.Zmin) * inv_log_q + 1e-12));

    const double probability_cost = max((double)0.0, -label.P);
    long long kP = 0;
    if (probability_cost > dpp.Pmin)
        kP = max(0LL, (long long)floor(
            log(probability_cost / dpp.Pmin) * inv_log_q + 1e-12));

    return {kW, kP};
}

void WernerAlgo_routing::bucket_by_ZP(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    map<pair<long long,long long>,ZLabel> buckets;
    for(const auto& L:cand){
        auto key = bucket_key(L);
        if(buckets.count(key)==0||L.B<buckets[key].B)
            buckets[key]=L;
    }
    vector<ZLabel> bucketed;
    for(const auto& L:buckets)
        bucketed.push_back(L.second);
    //pareto_prune_byZ(bucketed);
    // 排序鍵 = J 的 log 主體 Z²−P（eval_best_J 的 J=(α+B)·exp(Z²−P)，同 cell
    // 同 t 的 α、B 差異有限）。舊版按 Z 遞增排，Z 最小 = purify 最兇的 label
    // 永遠佔住 merge 的 per-k 前綴，導致「不 purify 也達標」的組合在中間層就
    // 被剪光，rounding 只剩 over-purify 的 shape（prob 崩、memory 爆）。
    sort(bucketed.begin(), bucketed.end(), [](const ZLabel& x, const ZLabel& y){
        return x.Z * x.Z - x.P < y.Z * y.Z - y.P;
    });
    cand.swap(bucketed);
}

Shape_vector WernerAlgo_routing::backtrack_shape(ZLabel leaf, vector<int>& out_purify_rounds){
    // a,b,k 已是節點 id（all-pairs DP），不再經過 path[] 映射
    int left_id=leaf.a,right_id=leaf.b;
    if(leaf.op==Op::LEAF){
        Shape_vector result;
        if (leaf.ent_l < 0 || leaf.ent_r < 0) return Shape_vector{};
        // 標準區間格式：[ent_l, ent_r]（跟原版 WernerAlgo 一致）
        result.push_back({left_id,  {{leaf.ent_l, leaf.ent_r}}});
        result.push_back({right_id, {{leaf.ent_l, leaf.ent_r}}});
        // LEAF: 一條 link，purify rounds = purify_type (即 tlen-1)
        out_purify_rounds = { leaf.purify_type };
        return result;
    }
    if(leaf.op==Op::CONT){
        assert(leaf.parent_id>=0&&leaf.parent_id<DP_table[leaf.t-1][leaf.a][leaf.b].size());
        ZLabel pre_label=DP_table[leaf.t-1][leaf.a][leaf.b][leaf.parent_id];
        // CONT: idle 不改變 purify rounds，直接透傳
        Shape_vector last_time=backtrack_shape(pre_label,out_purify_rounds);
        if (last_time.empty()) return Shape_vector{};
        auto & prel=last_time.front().second[0],&prer=last_time.back().second[0];
        assert(last_time.front().first==leaf.a);
        assert(last_time.back().first==leaf.b);
        assert(prel.second==leaf.t-1);
        assert(prer.second==leaf.t-1);
        prel.second++;
        prer.second++;
        return last_time;
    }
    if(leaf.op==Op::MERGE){
        Shape_vector left_result,right_result,result;
        assert(leaf.k>=0);
        int k_id=leaf.k;
        // MERGE: 左右各自遞迴取 purify rounds，再串接
        vector<int> left_rounds, right_rounds;
        ZLabel left_leaf=DP_table[leaf.t-1][leaf.a][leaf.k][leaf.left_id];
        left_result=backtrack_shape(left_leaf,left_rounds);
        ZLabel right_leaf=DP_table[leaf.t-1][leaf.k][leaf.b][leaf.right_id];
        right_result=backtrack_shape(right_leaf,right_rounds);
        if(DEBUG) {
            assert(left_result.front().first == leaf.a);
            assert(left_result.front().second[0].second == leaf.t - 1);
            assert(left_result.front().second.size() == 1);
            assert(left_result.back().first == k_id);
            assert(right_result.front().first == k_id);
            assert(right_result.back().first == leaf.b);
            assert(right_result.back().second[0].second == leaf.t - 1);
            assert(left_result.back().second.size() == 1);
        }

        for(int i = 0; i < (int)left_result.size(); i++) {
            result.push_back(left_result[i]);
        }
        result.back().second.push_back(right_result.front().second.front());
        for(int i = 1; i < (int)right_result.size(); i++) {
            result.push_back(right_result[i]);
        }

        result.front().second[0].second++;
        result.back().second[0].second++;
        // 串接左右的 purify rounds
        out_purify_rounds = left_rounds;
        out_purify_rounds.insert(out_purify_rounds.end(), right_rounds.begin(), right_rounds.end());
        return result;
    }
    out_purify_rounds.clear();
    return Shape_vector{};
}
bool WernerAlgo_routing::has_duplicate_nodes(const Shape_vector& sh){
    set<int> seen;
    for(const auto& P : sh){
        if(!seen.insert(P.first).second) return true;
    }
    return false;
}
int WernerAlgo_routing::split_dis(int s,int d,const WernerAlgo_routing::ZLabel& L){
    if(L.op!=WernerAlgo_routing::Op::MERGE||L.k<0) return 1000000000;
    // swap 樹平衡度：k 到兩端 hop 距離的差（原版是 path index 中點，
    // all-pairs 版 s/d/k 是節點 id，要用 hop 距離才有意義）
    return abs(hop_dist[s][L.k]-hop_dist[L.k][d]);
}
pair<double,WernerAlgo_routing::ZLabel> WernerAlgo_routing::eval_best_J(int s, int d, int t, double alp){
    double bestJ=1e18;
    int bestdis=1000000000;
    int flag=0;
    ZLabel tmp={};
    for(const auto& L:DP_table[t][s][d]){
        double J=(alp+L.B)*exp(L.Z*L.Z-L.P);
        int dis=split_dis(s,d,L);
        if(J+EPS<bestJ||(fabs(J-bestJ)<=EPS&&dis<bestdis)){
            bestJ=J;
            tmp=L;
            bestdis=dis;
            flag=1;
        }
    }
    if(flag) return {bestJ,tmp};
    else return {INF,tmp};
}

void WernerAlgo_routing::run() {
    int round = 1;
    while (round-- && !requests.empty()) {
        variable_initialize();
        //cerr << "\033[1;31m"<< "[WernerAlgo's parameter] : "<< dpp.Zmin<<" "<<dpp.eps_bucket<<" "<<dpp.eta<< "\033[0m"<< endl;
        int it=0;
        double eps=1e-4;
        const int REUSE = 20;  // 每次 oracle 找到 shape 後重複灌 REUSE 次流量
        while (obj+eps < 1.0) {
            it++;
            Shape_vector shape=separation_oracle();
            if (shape.empty()) break;

            // 取得此 shape 對應的 purify rounds
            vector<int> cur_purify_rounds;
            if(shape_purify_map.count(shape))
                cur_purify_rounds = shape_purify_map[shape];

            // 計算 total_need_q（只算一次，REUSE 次共用）
            map<pair<int,int>, int> total_need_q;
            for(int i=0;i<(int)shape.size();i++){
                for(pair<int,int> usedtime:shape[i].second){
                    int start=usedtime.first,end=usedtime.second;
                    for(int t=start;t<=end;t++)
                        total_need_q[{shape[i].first, t}]++;
                }
            }
            for(int li=0; li<(int)shape.size()-1; li++){
                int rounds = (li < (int)cur_purify_rounds.size()) ? cur_purify_rounds[li] : 0;
                if(rounds <= 0) continue;
                int link_start = shape[li].second.back().first;
                int u = shape[li].first, v = shape[li+1].first;
                for(int ti=0; ti<=rounds+1; ti++){
                    int extra = (int)Purify_in_vt[rounds][rounds + 1 - ti] - 1;
                    if(extra <= 0) continue;
                    int t = link_start + ti;
                    if(t >= graph.get_time_limit()) continue;
                    total_need_q[{u, t}] += extra;
                    total_need_q[{v, t}] += extra;
                }
            }

            // 同一個 shape 灌 REUSE 次流量
            for(int reuse = 0; reuse < REUSE && obj+eps < 1.0; reuse++) {
                double q = 1.0;
                for(auto& P : total_need_q){
                    int node_id = P.first.first, t = P.first.second;
                    double theta = P.second;
                    double cap = graph.get_node_memory_at(node_id, t);
                    if(cap > 0) q = min(q, cap / theta);
                    else q = 0;
                }
                if(q<=1e-10) break;
                int req_idx=-1;
                for(int i=0;i<requests.size();i++){
                    int ln=shape.front().first,rn=shape.back().first;
                    if(requests[i]==make_pair(ln,rn)){
                        if(req_idx==-1||alpha[req_idx]>alpha[i]){
                            req_idx=i;
                        }
                    }
                }
                if(req_idx==-1) break;
                x[req_idx][shape]+=q;
                double ori=alpha[req_idx];
                alpha[req_idx]=alpha[req_idx]*(1+epsilon*q);
                obj+=(alpha[req_idx]-ori);
                for(auto& P : total_need_q){
                    int node_id = P.first.first, t = P.first.second;
                    double theta = P.second;
                    double original = beta[node_id][t];
                    if(graph.get_node_memory_at(node_id, t) == 0) {
                        beta[node_id][t] = INF;
                    } else {
                        beta[node_id][t] = beta[node_id][t] * (1 + epsilon * (q / (graph.get_node_memory_at(node_id, t) / theta)));
                    }
                    obj += (beta[node_id][t] - original) * graph.get_node_memory_at(node_id, t);
                }
            }
        }
        cerr << "[" << algorithm_name << "] LP done, " << it << " oracle calls" << endl;
        vector<pair<double, Shape_vector>> shapes;

        for(int i = 0; i < (int)requests.size(); i++) {
            for(auto P : x[i]) {
                shapes.push_back({P.second, P.first});
            }
        }

        sort(shapes.begin(), shapes.end(), [](pair<double, Shape_vector> left, pair<double, Shape_vector> right) {
            return left.first > right.first;
        });

        vector<bool> used(requests.size(), false);
        vector<int> finished;
        // per-request 最終結果（routing trace 用）：accepted / fail_fid / fail_mem / no_shape
        vector<string> req_outcome(requests.size(), "no_shape");

        // [診斷] 統計 rounding 結果
        int round_total = 0, round_purified = 0, round_nopurify = 0;
        int round_fail_fid = 0, round_fail_mem = 0, round_fail_purify_mem = 0;
        // 每條 link 的 purification 詳細資訊
        struct LinkPurifyDetail {
            int u, v;           // link 端點
            int rounds;         // purify 輪數 (0=不做)
            double raw_f;       // 原始 F_init
            double raw_w;       // 原始 Werner = (4F-1)/3
            double purified_w;  // purify 後的 Werner
            double purified_f;  // purify 後的 Fidelity = (3W+1)/4
        };
        // purification 前後 fidelity/prob 記錄
        struct PurifyLogEntry {
            int src, dst, hop;
            bool has_purify;
            double fid_before, werner_before, prob_before, fidprob_before;
            double fid_after, werner_after, prob_after, fidprob_after;
            vector<LinkPurifyDetail> link_details;
        };
        vector<PurifyLogEntry> purify_log_entries;
        // 所有「被接受」的 request 的 routing 路徑（不限有無 purify）
        struct AcceptedPathEntry {
            int req_idx, src, dst, hop;
            bool has_purify;
            double fid, prob;
            string path_str;  // 節點序列 + 每段 purify rounds
            string mem_str;   // 每個節點的 memory 區間（可看出 swap 樹）
        };
        vector<AcceptedPathEntry> accepted_path_entries;

        for(pair<double, Shape_vector> P : shapes) {
            // 用 purify_rounds 構造 Shape（若有的話）
            vector<int> pr;
            if(shape_purify_map.count(P.second))
                pr = shape_purify_map[P.second];
            Shape shape = pr.empty() ? Shape(P.second) : Shape(P.second, pr);
            bool has_purify = false;
            for (int r : pr) if (r > 0) has_purify = true;
            int request_index = -1;
            for(int i = 0; i < (int)requests.size(); i++) {
                if(used[i] == false && requests[i] == make_pair(shape.get_node_mem_range().front().first, shape.get_node_mem_range().back().first)) {
                    request_index = i;
                }
            }

            if(request_index == -1 || used[request_index]) continue;
            round_total++;

            // 檢查資源（check_valid 可能 throw，用 try-catch 保護）
            bool resource_ok = false;
            bool fid_check = false;
            try {
                fid_check = graph.check_resource(shape, true, true);
            } catch(const runtime_error&) {
                fid_check = false;  // shape 不合法
            }
            if(!fid_check) {
                try {
                    bool mem_only = graph.check_resource(shape, false, true);
                    if(!mem_only) { round_fail_mem++; req_outcome[request_index] = "fail_mem"; }
                    else { round_fail_fid++; req_outcome[request_index] = "fail_fid"; }
                } catch(const runtime_error&) {
                    round_fail_fid++;
                    req_outcome[request_index] = "fail_fid";
                }
            }
            if(fid_check) {
                resource_ok = true;
                // 計算含 purification 的每 (node, timeslot) 總 memory 需求
                Shape_vector sv_chk = shape.get_node_mem_range();
                vector<int> pr_chk = shape.get_link_purify_rounds();
                map<pair<int,int>, int> total_need; // (node, t) -> total amount

                // 基本需求（跟 check_resource 一樣）
                for(size_t i = 0; i < sv_chk.size(); ++i) {
                    int node = sv_chk[i].first;
                    for(auto& rng : sv_chk[i].second) {
                        for(int t = rng.first; t <= rng.second; ++t)
                            total_need[{node, t}]++;
                    }
                }
                // purification 額外需求（index 需反向對應 gen_leaf_label）
                for(size_t li = 0; li < sv_chk.size() - 1; ++li) {
                    int rounds = (li < pr_chk.size()) ? pr_chk[li] : 0;
                    if(rounds <= 0) continue;
                    int link_start = sv_chk[li].second.back().first;
                    int u = sv_chk[li].first, v = sv_chk[li+1].first;
                    for(int ti = 0; ti <= rounds + 1; ++ti) {
                        int extra = (int)Purify_in_vt[rounds][rounds + 1 - ti] - 1;
                        if(extra <= 0) continue;
                        int t = link_start + ti;
                        if(t >= graph.get_time_limit()) { resource_ok = false; break; }
                        total_need[{u, t}] += extra;
                        total_need[{v, t}] += extra;
                    }
                    if(!resource_ok) break;
                }
                // 檢查所有 (node, t) 是否有足夠 memory
                for(auto& [nt, amount] : total_need) {
                    if(!resource_ok) break;
                    if(graph.get_node_memory_at(nt.first, nt.second) < amount)
                        resource_ok = false;
                }
            }
            if(!resource_ok && fid_check) {
                round_fail_purify_mem++;  // fidelity ok 但 purify memory 不夠
                req_outcome[request_index] = "fail_mem";
            }
            if(resource_ok) {
                used[request_index] = true;
                req_outcome[request_index] = "accepted";
                if(has_purify) round_purified++; else round_nopurify++;
                graph.reserve_shape(shape, true);
                // 額外扣除 purification 多消耗的 memory（標準 Shape 已扣 1，這裡補扣剩餘）
                // index 需反向對應 gen_leaf_label
                {
                    Shape_vector sv_res = shape.get_node_mem_range();
                    vector<int> pr_res = shape.get_link_purify_rounds();
                    for(size_t li = 0; li < sv_res.size() - 1; ++li) {
                        int rounds = (li < pr_res.size()) ? pr_res[li] : 0;
                        if(rounds <= 0) continue;
                        int link_start = sv_res[li].second.back().first;
                        int u = sv_res[li].first, v = sv_res[li+1].first;
                        for(int ti = 0; ti <= rounds + 1; ++ti) {
                            int extra = (int)Purify_in_vt[rounds][rounds + 1 - ti] - 1;
                            if(extra <= 0) continue;
                            int t = link_start + ti;
                            graph.reserve_node_memory_at(u, t, extra);
                            graph.reserve_node_memory_at(v, t, extra);
                        }
                    }
                }
                finished.push_back(request_index);

                // [新增] 記錄 routing 出來的路徑（node 序列 + 每段 purify rounds），
                // 同時印到 stderr 與收集到 accepted_path_entries 供寫檔
                {
                    Shape_vector sv_path = shape.get_node_mem_range();
                    vector<int> pr_path = shape.get_link_purify_rounds();
                    ostringstream path_ss, mem_ss;
                    for(size_t i = 0; i < sv_path.size(); ++i) {
                        path_ss << sv_path[i].first;
                        if(i + 1 < sv_path.size()) {
                            int r = (i < pr_path.size()) ? pr_path[i] : 0;
                            path_ss << " --(pur=" << r << ")--> ";
                        }
                        mem_ss << sv_path[i].first;
                        for(const auto& rng : sv_path[i].second)
                            mem_ss << "[" << rng.first << "," << rng.second << "]";
                        if(i + 1 < sv_path.size()) mem_ss << " ";
                    }
                    cerr << "\033[1;36m" << "[" << algorithm_name << " path] req#" << request_index
                         << " (" << sv_path.front().first << "->" << sv_path.back().first << "): "
                         << path_ss.str()
                         << " | hop=" << (sv_path.size() - 1) << "\033[0m" << endl;

                    double fid_acc = shape.get_fidelity(A, B, n, T, tao, graph.get_F_init(), true);
                    double prob_acc = has_purify ? graph.path_Pr_purify(shape) : graph.path_Pr(shape);
                    accepted_path_entries.push_back({
                        request_index,
                        sv_path.front().first, sv_path.back().first,
                        (int)sv_path.size() - 1,
                        has_purify, fid_acc, prob_acc,
                        path_ss.str(), mem_ss.str()
                    });
                    // routing trace CSV（per-request 對照用）
                    log_routing_trace(request_index, sv_path.front().first, sv_path.back().first,
                                      "accepted", sv_path, pr_path, fid_acc, prob_acc);
                }

                // [新增] 收集 purification 前後的 fidelity 與 prob 統計（僅記錄有 purify 的）
                if(has_purify) {
                    // 計算不開放 purify 的 fidelity 和 prob
                    Shape shape_no_pur(P.second);  // 不帶 purify rounds 的 shape
                    double fid_no_pur = shape_no_pur.get_fidelity(A, B, n, T, tao, graph.get_F_init(), false);
                    double prob_no_pur = graph.path_Pr(shape_no_pur);

                    // 計算開放 purify 後的 fidelity 和 prob
                    double fid_with_pur = shape.get_fidelity(A, B, n, T, tao, graph.get_F_init(), true);
                    double prob_with_pur = has_purify ? graph.path_Pr_purify(shape) : graph.path_Pr(shape);

                    Shape_vector sv = shape.get_node_mem_range();
                    int src = sv.front().first, dst = sv.back().first;
                    int hop = (int)sv.size() - 1;

                    // W = (4F - 1) / 3
                    double werner_no_pur = Purification::fidelity_to_werner(fid_no_pur);
                    double werner_with_pur = Purification::fidelity_to_werner(fid_with_pur);

                    // 收集每條 link 的 purify 詳細資訊
                    vector<LinkPurifyDetail> link_details;
                    vector<int> pur_rounds = shape.get_link_purify_rounds();
                    for (size_t li = 0; li < sv.size() - 1; ++li) {
                        int u = sv[li].first, v = sv[li+1].first;
                        int rounds = (li < pur_rounds.size()) ? pur_rounds[li] : 0;
                        double raw_f = graph.get_F_init(u, v);
                        double raw_w = Purification::fidelity_to_werner(raw_f);
                        double w_cur = Purification::pumping_werner(raw_w, rounds);
                        double purified_f = Purification::werner_to_fidelity(w_cur);
                        link_details.push_back({u, v, rounds, raw_f, raw_w, w_cur, purified_f});
                    }

                    purify_log_entries.push_back({
                        src, dst, hop, has_purify,
                        fid_no_pur, werner_no_pur, prob_no_pur, fid_no_pur * prob_no_pur,
                        fid_with_pur, werner_with_pur, prob_with_pur, fid_with_pur * prob_with_pur,
                        link_details
                    });
                }
            }
        }

        // [診斷] 印出 rounding 統計
        cerr << "\033[1;35m" << "[" << algorithm_name << " rounding] "
             << "candidates=" << round_total
             << " | accepted_purified=" << round_purified
             << " | accepted_nopurify=" << round_nopurify
             << " | fail_fidelity=" << round_fail_fid
             << " | fail_base_mem=" << round_fail_mem
             << " | fail_purify_extra_mem=" << round_fail_purify_mem
             << " | purify_map_size=" << shape_purify_map.size()
             << "\033[0m" << endl;

        // routing trace：所有未被接受的 request 也各記一行（NA 欄位）
        for(int i = 0; i < (int)requests.size(); i++) {
            if(used[i]) continue;
            log_routing_trace(i, requests[i].first, requests[i].second, req_outcome[i],
                              Shape_vector{}, vector<int>{}, 0, 0);
        }

        sort(finished.rbegin(), finished.rend());
        for(auto fin : finished) {
            requests.erase(requests.begin() + fin);
        }

        // [新增] 將 purification 前後比較統計寫入檔案 (append 模式, omp critical 保護)
        #pragma omp critical(zfa_routing_log_write)
        {
            string log_file_path = "../data/log/ZFA_routing_Purification_Stats.txt";
            ofstream log_file(log_file_path, ios::app);

            if (log_file.is_open()) {
                if (!experiment_label.empty()) {
                    log_file << "=== Experiment: " << experiment_label << " ===" << endl;
                }
                // 所有被接受的 routing 路徑（含無 purify 的）
                // mem 行：每個節點的 memory 佔用區間，interior 節點兩段區間的
                // 結束時間即該節點做 swap 的時刻，可據此還原 swap 樹
                log_file << "--- Accepted Paths (total: " << accepted_path_entries.size() << ") ---" << endl;
                for (auto& e : accepted_path_entries) {
                    log_file << "  req#" << e.req_idx
                             << " SD=(" << e.src << "," << e.dst << ")"
                             << " hop=" << e.hop
                             << " purified=" << (e.has_purify ? "YES" : "NO")
                             << " fid=" << e.fid
                             << " prob=" << e.prob << endl;
                    log_file << "    path: " << e.path_str << endl;
                    log_file << "    mem : " << e.mem_str << endl;
                }
                // 以下僅列「有做 purification」的 request 的前後比較
                log_file << "--- Purification Before/After (purified only, total: " << purify_log_entries.size() << ") ---" << endl;
                for (auto& e : purify_log_entries) {
                    log_file << "  SD=(" << e.src << "," << e.dst << ") hop=" << e.hop
                             << " purified=" << (e.has_purify ? "YES" : "NO") << endl;
                    log_file << "    [Before Purify] fidelity=" << e.fid_before
                             << "  werner=" << e.werner_before
                             << "  prob=" << e.prob_before
                             << "  fid*prob=" << e.fidprob_before << endl;
                    log_file << "    [After  Purify] fidelity=" << e.fid_after
                             << "  werner=" << e.werner_after
                             << "  prob=" << e.prob_after
                             << "  fid*prob=" << e.fidprob_after << endl;
                    // 每條 link 的 purify 詳細資訊
                    for (auto& lk : e.link_details) {
                        log_file << "      link(" << lk.u << "->" << lk.v << ")"
                                 << " purify_rounds=" << lk.rounds
                                 << "  F: " << lk.raw_f << " -> " << lk.purified_f
                                 << "  W: " << lk.raw_w << " -> " << lk.purified_w
                                 << endl;
                    }
                }
                log_file << "-----------------" << endl;
                log_file.close();
            } else {
                cerr << "[Warning] Unable to open log file: " << log_file_path << endl;
            }
        }
    }
    update_res();
    cerr << "[" << algorithm_name << "] end" << endl;
}
