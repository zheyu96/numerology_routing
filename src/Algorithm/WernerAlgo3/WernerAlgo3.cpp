#include "WernerAlgo3.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

WernerAlgo3::WernerAlgo3(const Graph& graph,const vector<pair<int,int>>& requests,const map<SDpair, vector<Path>>& paths): AlgorithmBase(graph, requests, paths)
{
    algorithm_name = "ZFA_UB";
}

// === 以下 LP 階段的函式和 WernerAlgo2 完全相同 ===

void WernerAlgo3::variable_initialize() {
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

Shape_vector WernerAlgo3::separation_oracle(){
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

WernerAlgo3::ZLabel WernerAlgo3::gen_leaf_label(int s,int e,int st,int tlen,int path_a,int path_b) {
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

void WernerAlgo3::run_dp_in_t(const Path& path, const DPParam& dpp,int t) {
    const int T = graph.get_time_limit();
    const int n = (int)path.size();
    const size_t MAX_CANDIDATES_PER_CELL = 1000;
    const size_t MAX_LABELS_PER_CELL = 100;
    auto trim_cell = [&](vector<ZLabel>& labels) {
        if(labels.size() <= MAX_LABELS_PER_CELL) return;
        nth_element(labels.begin(), labels.begin() + MAX_LABELS_PER_CELL, labels.end(),
            [](const ZLabel& x, const ZLabel& y) {
                if(x.B != y.B) return x.B < y.B;
                return x.Z < y.Z;
            });
        labels.resize(MAX_LABELS_PER_CELL);
        bucket_by_ZP(labels);
        if(labels.size() > MAX_LABELS_PER_CELL) {
            sort(labels.begin(), labels.end(), [](const ZLabel& x, const ZLabel& y) {
                if(x.B != y.B) return x.B < y.B;
                return x.Z < y.Z;
            });
            labels.resize(MAX_LABELS_PER_CELL);
        }
    };
    for(int a=0;a<n-1;a++)
        for(int b=a+1;b<n;b++){
            int s=path[a],e=path[b];
            vector<ZLabel> cand;
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
            const auto& pre=DP_table[t-1][a][b];
            for(int p_id=0;p_id<pre.size() && cand.size() < MAX_CANDIDATES_PER_CELL;p_id++){
                double Zp=pre[p_id].Z+dpp.eta;
                if(Zp<=dpp.Zhat){
                    double Bp=pre[p_id].B+beta[s][t]+beta[e][t];
                    double Pp=pre[p_id].P;
                    ZLabel L(Bp,Zp,Pp,Op::CONT,-1,a,b,t,-1,p_id);
                    cand.push_back(L);
                }
            }
            for(int k=a+1;k<b;k++){
                const auto& L1=DP_table[t-1][a][k];
                const auto& L2=DP_table[t-1][k][b];
                if(L1.size()==0||L2.size()==0) continue;
                for(int lid=0;lid<L1.size() && cand.size() < MAX_CANDIDATES_PER_CELL;lid++)
                    for(int rid=0;rid<L2.size() && cand.size() < MAX_CANDIDATES_PER_CELL;rid++){
                        const auto& left_seg=L1[lid];
                        const auto& right_seg=L2[rid];
                        double Zp=sqrt((left_seg.Z+dpp.eta)*(left_seg.Z+dpp.eta)+
                                        (right_seg.Z+dpp.eta)*(right_seg.Z+dpp.eta));
                        double swap_prob=log(graph.get_node_swap_prob(path[k]));
                        double Pp=left_seg.P+right_seg.P+swap_prob;
                        if(Zp<=dpp.Zhat){
                            double Bp=left_seg.B+right_seg.B+beta[s][t]+beta[e][t];
                            ZLabel L(Bp,Zp,Pp,Op::MERGE,-1,a,b,t,k,-1,lid,rid);
                            cand.push_back(L);
                        }
                    }
            }
            vector<ZLabel> non_leaf;
            for (auto& L : cand) {
                if (L.op != Op::LEAF)
                    non_leaf.push_back(L);
            }
            bucket_by_ZP(non_leaf);
            trim_cell(non_leaf);
            cand.erase(
                remove_if(cand.begin(), cand.end(),
                        [](const ZLabel& L){ return L.op != Op::LEAF; }),
                cand.end());
            cand.insert(cand.end(), non_leaf.begin(), non_leaf.end());
            trim_cell(cand);
            DP_table[t][a][b] = cand;
        }
}

void WernerAlgo3::pareto_prune_byZ(vector<ZLabel>& cand) {
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

void WernerAlgo3::bucket_by_ZP(vector<ZLabel>& cand) {
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
    sort(bucketed.begin(), bucketed.end(), [](const ZLabel& x, const ZLabel& y){
        return x.Z < y.Z;
    });
    cand.swap(bucketed);
}

Shape_vector WernerAlgo3::backtrack_shape(ZLabel leaf,const vector<int>& path, vector<int>& out_purify_rounds){
    int left_id=path[leaf.a],right_id=path[leaf.b];
    if(leaf.op==Op::LEAF){
        Shape_vector result;
        if (leaf.ent_time.size() < 2) return Shape_vector{};
        assert(leaf.ent_time.size() == 2);
        result.push_back({left_id,  {{leaf.ent_time[0], leaf.ent_time[1]}}});
        result.push_back({right_id, {{leaf.ent_time[0], leaf.ent_time[1]}}});
        out_purify_rounds = { leaf.purify_type };
        return result;
    }
    if(leaf.op==Op::CONT){
        assert(leaf.parent_id>=0&&leaf.parent_id<DP_table[leaf.t-1][leaf.a][leaf.b].size());
        ZLabel pre_label=DP_table[leaf.t-1][leaf.a][leaf.b][leaf.parent_id];
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
        out_purify_rounds = left_rounds;
        out_purify_rounds.insert(out_purify_rounds.end(), right_rounds.begin(), right_rounds.end());
        return result;
    }
    out_purify_rounds.clear();
    return Shape_vector{};
}

int WernerAlgo3::split_dis(int s,int d,WernerAlgo3::ZLabel& L){
    if(L.op!=WernerAlgo3::Op::MERGE||L.k<0) return 1000000000;
    int mid=(s+d)/2;
    return abs(mid-L.k);
}

pair<double,WernerAlgo3::ZLabel> WernerAlgo3::eval_best_J(int s, int d, int t, double alp){
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

// === run()：LP 階段同 ZFA2，但不做 greedy rounding，改算 LP upper bound ===

void WernerAlgo3::run() {
    int round = 1;
    while (round-- && !requests.empty()) {
        variable_initialize();
        int it=0;
        double eps=1e-4;
        const int REUSE = 20;
        cerr << "[" << algorithm_name << "] LP phase start, requests=" << requests.size() << endl;
        while (obj+eps < 1.0) {
            it++;
            Shape_vector shape=separation_oracle();
            if (shape.empty()) { cerr << "[" << algorithm_name << "] iter " << it << ": oracle empty, stop" << endl; break; }
            vector<int> cur_purify_rounds;
            if(shape_purify_map.count(shape))
                cur_purify_rounds = shape_purify_map[shape];

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
                if (it % 10 == 1 && reuse == 0) {
                    int src=shape.front().first, dst=shape.back().first;
                    cerr << "[" << algorithm_name << "] oracle#" << it
                         << " | obj=" << (double)obj
                         << " | q=" << (double)q
                         << " | req=(" << src << "," << dst << ")"
                         << " | hops=" << (shape.size()-1)
                         << endl;
                }
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

        // === LP Upper Bound（論文 Eq. 14a-14d）===
        // Eq.14a: max Σ_{i,p,m} Pr(i,p,m) · w(i,p,m) · x_m^ip
        // Eq.14b: Σ_{p,m} x_m^ip ≤ 1, ∀i ∈ I
        // Garg-Konemann 累積的 x 可能違反 14b（Σ x > 1），需 per-request 正規化
        cerr << "[" << algorithm_name << "] LP done, " << it << " iters, obj=" << (double)obj << endl;
        int ub_cnt = 0;
        double max_xsum = 0;
        for(int i = 0; i < (int)requests.size(); i++) {
            if(x[i].empty()) continue;
            // 計算此 request 的 Σ x_m^ip
            double xsum = 0;
            for(auto& P : x[i]) xsum += P.second;
            if(xsum < 1e-12) continue;
            max_xsum = max(max_xsum, (double)xsum);
            // Eq.14b 正規化: scale down 使 Σ x ≤ 1
            double scale = (xsum > 1.0L) ? (1.0L / xsum) : 1.0L;

            for(auto& P : x[i]) {
                double xval = P.second * scale;
                if(xval < 1e-12) continue;
                vector<int> pr;
                if(shape_purify_map.count(P.first))
                    pr = shape_purify_map[P.first];
                Shape shape = pr.empty() ? Shape(P.first) : Shape(P.first, pr);
                double fidelity;
                try {
                    fidelity = shape.get_fidelity(A, B, n, T, tao, graph.get_F_init(), true);
                } catch(const runtime_error&) {
                    continue;  // shape timing 不合法，跳過
                }
                if(fidelity + EPS < graph.get_fidelity_threshold()) continue;
                double Pr = graph.path_Pr_purify(shape);
                double w = (4.0L * fidelity - 1.0L) / 3.0L;
                res["fidelity_gain"] += xval * Pr * w;
                res["succ_request_cnt"] += xval * Pr;
                ub_cnt++;
            }
        }
        // 將 LP upper bound 放大 1.1 倍，作為更寬鬆的上界
        res["fidelity_gain"] *= 1.1;
        res["succ_request_cnt"] *= 1.1;
        res["actual_req_cnt"] = (double)ub_cnt;

        cerr << "[" << algorithm_name << "] UB(Eq.14a)*1.1: fid_gain=" << (double)res["fidelity_gain"]
             << " succ_req=" << (double)res["succ_request_cnt"]
             << " actual_req=" << (double)res["actual_req_cnt"]
             << " shapes=" << ub_cnt
             << " max_xsum=" << max_xsum << endl;
    }
    cerr << "[" << algorithm_name << "] end" << endl;
}
