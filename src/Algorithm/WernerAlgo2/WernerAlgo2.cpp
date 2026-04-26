#include "WernerAlgo2.h"
#include <fstream> 
#include <iostream>
#include <cmath>

using namespace std;

WernerAlgo2::WernerAlgo2(const Graph& graph,const vector<pair<int,int>>& requests,const map<SDpair, vector<Path>>& paths): AlgorithmBase(graph, requests, paths)
{
    algorithm_name = "ZFA2";
}

void WernerAlgo2::variable_initialize() {
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
    dpp.eps_bucket = graph.get_bucket_eps();
    double F_th=graph.get_fidelity_threshold();
    double w_th=(4.0*F_th-1.0)/3.0;
    dpp.Zhat = sqrt(-log(w_th))+1e-9;
    dpp.Zmin = graph.get_Zmin();
    dpp.T    = time_limit-1;
    dpp.tau_max=min(time_limit-1,5);
    dpp.eta  = graph.get_tao()/graph.get_time_limit();
    dpp.deltaP=graph.get_delta_P();
    beta.assign(V, vector<double>(T, INF));

    for (int v = 0; v < V; ++v) {
        for (int t = 0; t < T; ++t) {
            int cap = graph.get_node_memory_at(v, t);
            beta[v][t] = (cap == 0) ? INF : (delta / cap);
        }
    }

    oracle_cache.clear();
    oracle_cache.resize(requests.size());
    for (int i = 0; i < (int)requests.size(); i++)
        oracle_cache[i].resize(get_paths(requests[i].first, requests[i].second).size());
    dirty_nodes.clear();
    dirty_alpha_idxs.clear();
}

Shape_vector WernerAlgo2::separation_oracle(){
    double most_violate=1e9;
    Shape_vector todo_shape;
    vector<int> best_purify_rounds;

    for(int i=0;i<(int)requests.size();i++){
        int src=requests[i].first,dst=requests[i].second;
        const vector<Path>& cur_paths=get_paths(src,dst);

        for(int p=0;p<(int)cur_paths.size();p++){
            auto& cache = oracle_cache[i][p];

            if(cache.valid && !dirty_alpha_idxs.count(i)){
                bool path_dirty = false;
                for(int v : cur_paths[p]){
                    if(dirty_nodes.count(v)){ path_dirty = true; break; }
                }
                if(!path_dirty){
                    if(cache.best_score < most_violate){
                        most_violate = cache.best_score;
                        todo_shape = cache.shape;
                        best_purify_rounds = cache.purify_rounds;
                    }
                    continue;
                }
            }

            // Full recompute
            int T=dpp.T+5;
            int n=cur_paths[p].size()+5;
            DP_table.clear();
            DP_table.resize(T);
            for(int ii=0;ii<(int)DP_table.size();ii++){
                DP_table[ii].resize(n);
                for(int j=0;j<(int)DP_table[ii].size();j++)
                    DP_table[ii][j].resize(n);
            }

            double local_best_J = 1e18;
            ZLabel local_best_label;

            for(int t=1;t<=dpp.T;t++){
                run_dp_in_t(cur_paths[p],dpp,t);
                auto cur_val=eval_best_J(0,cur_paths[p].size()-1,t,alpha[i]);
                if(cur_val.first < local_best_J){
                    local_best_J = cur_val.first;
                    local_best_label = cur_val.second;
                }
            }

            if(local_best_J < 1e18){
                vector<int> cur_rounds;
                cache.shape = backtrack_shape(local_best_label, cur_paths[p], cur_rounds);
                cache.purify_rounds = cur_rounds;
                cache.best_score = local_best_J;
                cache.valid = true;

                if(local_best_J < most_violate){
                    most_violate = local_best_J;
                    todo_shape = cache.shape;
                    best_purify_rounds = cur_rounds;
                }
            }
        }
    }

    dirty_nodes.clear();
    dirty_alpha_idxs.clear();

    if(!todo_shape.empty()){
        auto it = shape_purify_map.find(todo_shape);
        if(it == shape_purify_map.end()){
            shape_purify_map[todo_shape] = best_purify_rounds;
        } else {
            bool new_has = false, old_has = false;
            for(int r : best_purify_rounds) if(r > 0) new_has = true;
            for(int r : it->second) if(r > 0) old_has = true;
            if(new_has && !old_has)
                it->second = best_purify_rounds;
        }
    }
    return todo_shape;
}
WernerAlgo2::ZLabel WernerAlgo2::gen_leaf_label(int s,int e,int st,int tlen,int path_a,int path_b) {
    double Bleaf=0.0;
    if(st-tlen<0) return ZLabel();
    for(int i=0;i<=tlen;i++){
        double bt=beta[s][st-i]+beta[e][st-i];
        Bleaf+=bt*Purify_in_vt[tlen-1][i];
    }
    // Werner-equivalent of paper Eq.(4)-(5) under F=(3w+1)/4:
    //   P_p   = (w1*w2 + 1) / 2
    //   w_new = (4*w1*w2 + w1 + w2) / (3*w1*w2 + 3)
    double w_ini=graph.get_link_werner(s,e);
    double w_cur=w_ini,p_cur=graph.get_entangle_succ_prob(s,e);
    for(int i=1;i<=tlen-1;i++){
        p_cur*=(w_cur*w_ini+1.0L)/2.0L;
        w_cur=(4.0L*w_cur*w_ini+w_cur+w_ini)/(3.0L*w_cur*w_ini+3.0L);
        p_cur*=graph.get_entangle_succ_prob(s,e);
    }
    double Zleaf=sqrt(-log(w_cur));
    double Pleaf=log(p_cur);
    if(Zleaf>dpp.Zhat) return ZLabel();
    return ZLabel(Bleaf,Zleaf,Pleaf,Op::LEAF,tlen-1,path_a,path_b,st,-1);
} 
void WernerAlgo2::run_dp_in_t(const Path& path, const DPParam& dpp,int t) {
    const int T = graph.get_time_limit();
    const int n = (int)path.size();

    // -------- t = 1..T-1 外圈時間迴圈 --------
    for(int a=0;a<n-1;a++)
        for(int b=a+1;b<n;b++){
            int s=path[a],e=path[b];
            vector<ZLabel> cand;
            //leaf
            if(a+1==b){
                for(int i=0;i<=purify_time;i++){
                    if(t-i-1<=0) continue;
                    ZLabel L=gen_leaf_label(s,e,t,i+1,a,b);
                    if(L.Z<=dpp.Zhat){
                        L.ent_time={t-i-1,t};
                        cand.push_back(L);
                    }
                }
            }
            //continue
            auto pre=DP_table[t-1][a][b];
            for(int p_id=0;p_id<pre.size();p_id++){
                double Zp=pre[p_id].Z+dpp.eta;
                if(Zp<=dpp.Zhat){
                    double Bp=pre[p_id].B+beta[s][t]+beta[e][t];
                    double Pp=pre[p_id].P;
                    ZLabel L(Bp,Zp,Pp,Op::CONT,-1,a,b,t,-1,p_id);
                    cand.push_back(L);
                }
            }
            //merge
            for(int k=a+1;k<b;k++){
                auto L1=DP_table[t-1][a][k],L2=DP_table[t-1][k][b];
                if(L1.size()==0||L2.size()==0) continue;
                for(int lid=0;lid<L1.size();lid++)
                    for(int rid=0;rid<L2.size();rid++){
                        auto left_seg=L1[lid],right_seg=L2[rid];
                        double Zp=sqrt((left_seg.Z+dpp.eta)*(left_seg.Z+dpp.eta)+
                                        (right_seg.Z+dpp.eta)*(right_seg.Z+dpp.eta));
                        double swap_prob=log(graph.get_node_swap_prob(path[k]));
                        double Pp=left_seg.P+right_seg.P+swap_prob;
                        if(Zp<=dpp.Zhat){
                            double Bp=left_seg.B+right_seg.B+beta[s][t]+beta[e][t];
                            ZLabel L(Bp,Zp,Pp,Op::MERGE,-1,a,b,t,k,-1,lid,rid);
                            cand.push_back(L);
                            // debug: 只印第一次 merge 的 Z 細節
                            static bool merge_printed = false;
                            if(!merge_printed) {
                                cerr << "[ZFA2:merge] PASS Zp=" << (double)Zp << " Zhat=" << (double)dpp.Zhat
                                     << " leftZ=" << (double)left_seg.Z << " rightZ=" << (double)right_seg.Z
                                     << " leftOp=" << (int)left_seg.op << " rightOp=" << (int)right_seg.op
                                     << " leftPurify=" << left_seg.purify_type << " rightPurify=" << right_seg.purify_type
                                     << " t=" << t << endl;
                                merge_printed = true;
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
            cand.erase(
                remove_if(cand.begin(), cand.end(),
                        [](const ZLabel& L){ return L.op != Op::LEAF; }),
                cand.end());
            cand.insert(cand.end(), non_leaf.begin(), non_leaf.end());
            DP_table[t][a][b] = cand;
        }
}

void WernerAlgo2::pareto_prune_byZ(vector<ZLabel>& cand) {
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

void WernerAlgo2::bucket_by_ZP(vector<ZLabel>& cand) {
    if (cand.empty()) return;
    double q=1+dpp.eps_bucket;
    double invLogQ=1.0/log(q);
    double deltaP=dpp.deltaP;
    map<pair<double,double>,ZLabel> buckets;
    for(auto L:cand){
        double kW;
        if(L.Z<=dpp.Zmin) kW=0.0;
        else{
            kW=floor(log(L.Z/dpp.Zmin)*invLogQ+1e-12);
            if(kW<0) kW=0.0;
        } 
        double kP=floor(-L.P/deltaP);
        auto key=make_pair(kW,kP);
        if(buckets.count(key)==0||L.B<buckets[key].B)
            buckets[key]=L;
    }
    vector<ZLabel> bucketed;
    for(auto L:buckets)
        bucketed.push_back(L.second);
    //pareto_prune_byZ(bucketed);
    sort(bucketed.begin(), bucketed.end(), [](const ZLabel& x, const ZLabel& y){
        return x.Z < y.Z;
    });
    cand.swap(bucketed);
}

Shape_vector WernerAlgo2::backtrack_shape(ZLabel leaf,const vector<int>& path, vector<int>& out_purify_rounds){
    int left_id=path[leaf.a],right_id=path[leaf.b];
    if(leaf.op==Op::LEAF){
        Shape_vector result;
        if (leaf.ent_time.size() < 2) return Shape_vector{};
        assert(leaf.ent_time.size() == 2);
        // 標準區間格式：[ent_time[0], ent_time[1]]（跟原版 WernerAlgo 一致）
        result.push_back({left_id,  {{leaf.ent_time[0], leaf.ent_time[1]}}});
        result.push_back({right_id, {{leaf.ent_time[0], leaf.ent_time[1]}}});
        // LEAF: 一條 link，purify rounds = purify_type (即 tlen-1)
        out_purify_rounds = { leaf.purify_type };
        return result;
    }
    if(leaf.op==Op::CONT){
        assert(leaf.parent_id>=0&&leaf.parent_id<DP_table[leaf.t-1][leaf.a][leaf.b].size());
        ZLabel pre_label=DP_table[leaf.t-1][leaf.a][leaf.b][leaf.parent_id];
        // CONT: idle 不改變 purify rounds，直接透傳
        Shape_vector last_time=backtrack_shape(pre_label,path,out_purify_rounds);
        if (last_time.empty()) return Shape_vector{};
        auto & prel=last_time.front().second[0],&prer=last_time.back().second[0];
        assert(last_time.front().first==path[leaf.a]);
        assert(last_time.back().first==path[leaf.b]);
        assert(prel.second==leaf.t-1);
        assert(prer.second==leaf.t-1);
        prel.second++;
        prer.second++;
        return last_time;
    }
    if(leaf.op==Op::MERGE){
        Shape_vector left_result,right_result,result;
        assert(leaf.k>=0);
        int k_id=path[leaf.k];
        // MERGE: 左右各自遞迴取 purify rounds，再串接
        vector<int> left_rounds, right_rounds;
        ZLabel left_leaf=DP_table[leaf.t-1][leaf.a][leaf.k][leaf.left_id];
        left_result=backtrack_shape(left_leaf,path,left_rounds);
        ZLabel right_leaf=DP_table[leaf.t-1][leaf.k][leaf.b][leaf.right_id];
        right_result=backtrack_shape(right_leaf,path,right_rounds);
        if(DEBUG) {
            assert(left_result.front().first == path[leaf.a]);
            assert(left_result.front().second[0].second == leaf.t - 1);
            assert(left_result.front().second.size() == 1);
            assert(left_result.back().first == k_id);
            assert(right_result.front().first == k_id);
            assert(right_result.back().first == path[leaf.b]);
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
int WernerAlgo2::split_dis(int s,int d,WernerAlgo2::ZLabel& L){
    if(L.op!=WernerAlgo2::Op::MERGE||L.k<0) return 1000000000;
    int mid=(s+d)/2;
    return abs(mid-L.k);
}
pair<double,WernerAlgo2::ZLabel> WernerAlgo2::eval_best_J(int s, int d, int t, double alp){
    double bestJ=1e18;
    int bestdis=1000000000;
    int flag=0;
    ZLabel tmp={};
    for(auto L:DP_table[t][s][d]){
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

void WernerAlgo2::run() {
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
                dirty_alpha_idxs.insert(req_idx);
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

            // Mark dirty nodes for incremental oracle
            for(int i=0;i<(int)shape.size();i++)
                dirty_nodes.insert(shape[i].first);
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
                    if(!mem_only) round_fail_mem++;
                    else round_fail_fid++;
                } catch(const runtime_error&) {
                    round_fail_fid++;
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
            if(!resource_ok && fid_check) round_fail_purify_mem++;  // fidelity ok 但 purify memory 不夠
            if(resource_ok) {
                used[request_index] = true;
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
                    double werner_no_pur = (4.0 * fid_no_pur - 1.0) / 3.0;
                    double werner_with_pur = (4.0 * fid_with_pur - 1.0) / 3.0;

                    // 收集每條 link 的 purify 詳細資訊
                    vector<LinkPurifyDetail> link_details;
                    vector<int> pur_rounds = shape.get_link_purify_rounds();
                    for (size_t li = 0; li < sv.size() - 1; ++li) {
                        int u = sv[li].first, v = sv[li+1].first;
                        int rounds = (li < pur_rounds.size()) ? pur_rounds[li] : 0;
                        double raw_f = graph.get_F_init(u, v);
                        double raw_w = (4.0 * raw_f - 1.0) / 3.0;
                        double w_cur = raw_w;
                        if (rounds > 0) {
                            // Werner-equivalent of paper Eq.(4):
                            //   w_new = (4*w1*w2 + w1 + w2) / (3*w1*w2 + 3)
                            for (int r = 0; r < rounds; r++) {
                                w_cur = (4.0L * w_cur * raw_w + w_cur + raw_w)
                                      / (3.0L * w_cur * raw_w + 3.0L);
                            }
                        }
                        double purified_f = (3.0 * w_cur + 1.0) / 4.0;
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

        sort(finished.rbegin(), finished.rend());
        for(auto fin : finished) {
            requests.erase(requests.begin() + fin);
        }

        // [新增] 將 purification 前後比較統計寫入檔案 (append 模式, omp critical 保護)
        #pragma omp critical(zfa2_log_write)
        {
            string log_file_path = "../data/log/ZFA2_Purification_Stats.txt";
            ofstream log_file(log_file_path, ios::app);

            if (log_file.is_open()) {
                if (!experiment_label.empty()) {
                    log_file << "=== Experiment: " << experiment_label << " ===" << endl;
                }
                log_file << "--- Accepted Requests (total: " << purify_log_entries.size() << ") ---" << endl;
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
