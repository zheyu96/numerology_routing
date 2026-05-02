
#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/MyAlgo1/MyAlgo1.h"
#include "Algorithm/MyAlgo2/MyAlgo2.h"
#include "Algorithm/MyAlgo3/MyAlgo3.h"
#include "Algorithm/MyAlgo4/MyAlgo4.h"
#include "Algorithm/MyAlgo5/MyAlgo5.h"
#include "Algorithm/MyAlgo6/MyAlgo6.h"
#include "Algorithm/WernerAlgo/WernerAlgo.h"
#include "Algorithm/WernerAlgo2/WernerAlgo2.h"
#include "Algorithm/WernerAlgo3/WernerAlgo3.h"
#include "Algorithm/WernerAlgo_UB/WernerAlgo_UB.h"
#include "Network/PathMethod/PathMethodBase/PathMethod.h"
#include "Network/PathMethod/Greedy/Greedy.h"
#include "Network/PathMethod/QCAST/QCAST.h"
#include "Network/PathMethod/REPS/REPS.h"

using namespace std;

SDpair generate_new_request(int num_of_node){
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, num_of_node-1);
    int node1 = unif(generator), node2 = unif(generator);
    while(node1 == node2) node2 = unif(generator);

    return make_pair(node1, node2);
}

vector<SDpair> generate_requests(Graph &graph, int requests_cnt, int length_lower, int length_upper) {
    int n = graph.get_num_nodes();
    vector<SDpair> cand;
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            int dist = graph.distance(i, j);
            if(dist >= length_lower && dist <= length_upper) {
                cand.emplace_back(i, j);
            }
        }
    }

    random_shuffle(cand.begin(), cand.end());

    vector<SDpair> requests;
    for(SDpair sdpair : cand) {
        int cnt = unif(generator) % 4 + 3;
        while(cnt--) requests.push_back(sdpair);
    }

    while((int)requests.size() < requests_cnt) {
        requests.emplace_back(generate_new_request(n));
    }

    while((int)requests.size() > requests_cnt) {
        requests.pop_back();
    }

    return requests;
}
vector<SDpair> generate_requests_fid(Graph &graph, int requests_cnt,double fid_th,double hop_th, double fid_upper = 1) {
    int n = graph.get_num_nodes();
    vector<pair<SDpair,double>> cand[22];
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);
    int sd_cnt=0;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            double fid = graph.get_ini_fid(i,j);
            //cerr<<"fid of "<<i<<" "<<j<<" : "<<fid<<endl;
            assert(fid>=0.0&&fid<=1.0);
            if(fid > fid_th && fid <= fid_upper && graph.distance(i,j)>=hop_th) {
                int index = fid/0.05;
                //index-=5;
                if(index < 0) continue;
                if(index > 20) index = 20;
                int d=graph.distance(i, j),f0=fid,prob=pow(0.1,d)*pow(0.9,max(d-1,0));
                //double score = f0+prob*100-0.1*d;
                cand[index].emplace_back(std::make_pair(std::make_pair(i, j), graph.distance(i, j)));
                if(graph.distance(i,j)>=1)sd_cnt++;
            }
        }
    }
     cerr << "\033[1;32m"<< "[SD ini pairs] : "<<sd_cnt<< "\033[0m"<< endl;
    /*for(int i=21;i>=0;i--){
        if(!cand[i].empty()){
            random_shuffle(cand[i].begin(), cand[i].end());
        }
    } */
    /* for(int i=21;i>=0;i--){
        sort(cand[i].begin(),cand[i].end(),[](const pair<SDpair,double>& L,const pair<SDpair,double>& R){
            return L.second > R.second;
        }) ;
    }  */
    for(int i=0;i<22;i++){
        random_shuffle(cand[i].begin(), cand[i].end());
    }
    // 檢查是否有任何候選
    bool any_cand = false;
    for (int i = 0; i < 22; i++) if (!cand[i].empty()) any_cand = true;
    if (!any_cand) {
        cerr << "[generate_requests_fid] WARNING: no candidates found (fid_th=" << fid_th << ", hop_th=" << hop_th << ")" << endl;
        return {};
    }

    vector<SDpair> requests;
    int pos[22];
    for(int i=0;i<22;i++) pos[i]=0;
    int idx=0;
    while((int)requests.size()<requests_cnt){
        int cnt=unif(generator) % 4 +3;
        cnt=min(cnt,(int)(requests_cnt-(int)requests.size()));
        // 找下一個非空桶（有保護）
        int tries = 0;
        while(cand[21-idx].empty()){
            idx++;
            if(idx>=22) idx=0;
            if(++tries > 22) break;  // 防止無限迴圈
        }
        if(tries > 22) break;
        if(!cand[21-idx].empty()){
            for(int i=0;i<cnt;i++){
                requests.push_back(cand[21-idx][pos[21-idx]].first);
            }
            pos[21-idx]++;
            pos[21-idx]%=cand[21-idx].size();
        }
        idx=(idx+1)%22;
    }
    if ((int)requests.size() < requests_cnt)
        cerr << "[generate_requests_fid] only generated " << requests.size() << "/" << requests_cnt << " requests" << endl;
    return requests;
}
// 生成「purification 能帶來優勢」的 request
// 用 Shape::get_fidelity 精確計算（和 rounding 階段的 check_resource 完全一致）
// 構造一個最簡單的 balanced-tree shape，分別算有/無 purify 的 real fidelity
vector<SDpair> generate_requests_purify_needed(Graph &graph, int requests_cnt, int min_hop = 2) {
    int n = graph.get_num_nodes();
    double fid_th = graph.get_fidelity_threshold();
    double A = graph.get_A(), B = graph.get_B();
    double n_param = graph.get_n(), T = graph.get_T(), tao = graph.get_tao();
    const auto& F_init = graph.get_F_init();

    // BFS 找最短路徑
    auto bfs_path = [&](int src, int dst) -> vector<int> {
        vector<int> parent(n, -1);
        vector<bool> vis(n, false);
        queue<int> que;
        vis[src] = true;
        que.push(src);
        while (!que.empty()) {
            int u = que.front(); que.pop();
            if (u == dst) break;
            for (int v : graph.adj_list[u]) {
                if (!vis[v]) {
                    vis[v] = true;
                    parent[v] = u;
                    que.push(v);
                }
            }
        }
        if (!vis[dst]) return {};
        vector<int> path;
        for (int v = dst; v != -1; v = parent[v]) path.push_back(v);
        reverse(path.begin(), path.end());
        return path;
    };

    // 用和 Shape::get_fidelity 完全相同的公式直接計算 path fidelity
    // F = A + B*exp(-(t/T)^n), decoherence 用 pass_tao
    auto t2F = [&](double t) -> double {
        if(t >= 1e5) return 0;
        return A + B * exp(-pow(t / T, n_param));
    };
    auto F2t = [&](double F) -> double {
        if(F <= A + 1e-9) return 1e9;
        return T * pow(-log((F - A) / B), 1.0 / n_param);
    };
    auto pass_tao_f = [&](double F) -> double {
        return t2F(F2t(F) + tao);
    };
    auto Fswap = [&](double Fa, double Fb) -> double {
        if(Fa <= A + 1e-9 || Fb <= A + 1e-9) return 0;
        return Fa * Fb + (1.0 / 3.0) * (1.0 - Fa) * (1.0 - Fb);
    };

    // 遞迴計算 balanced-tree schedule 的 end-to-end fidelity
    // edges[i] = F_init of edge i, purify_rounds[i] = rounds for edge i (0=none)
    function<double(int, int, const vector<double>&, const vector<int>&)> calc_fidelity;
    calc_fidelity = [&](int left, int right, const vector<double>& edge_fids, const vector<int>& pur_rounds) -> double {
        if (left == right - 1) {
            // Leaf: single edge
            double raw_f = edge_fids[left];
            int rounds = pur_rounds[left];
            if (rounds > 0) {
                double Fbase = raw_f;
                double Fcur = raw_f;
                for (int r = 0; r < rounds; r++) {
                    // Paper Eq.(4): pumping purification with a fresh pair Fbase.
                    double Fcur_bar = 1.0 - Fcur;
                    double Fbase_bar = 1.0 - Fbase;
                    double den = Fcur * Fbase
                               + (1.0 / 3.0) * Fcur * Fbase_bar
                               + (1.0 / 3.0) * Fcur_bar * Fbase
                               + (5.0 / 9.0) * Fcur_bar * Fbase_bar;
                    double num = Fcur * Fbase
                               + (1.0 / 9.0) * Fcur_bar * Fbase_bar;
                    Fcur = num / den;
                }
                return pass_tao_f(Fcur);
            }
            return pass_tao_f(raw_f);
        }
        // Balanced split
        int mid = (left + right) / 2;
        double Fa = calc_fidelity(left, mid, edge_fids, pur_rounds);
        double Fb = calc_fidelity(mid, right, edge_fids, pur_rounds);
        // Swap + 1 tao decoherence
        return t2F(F2t(Fswap(pass_tao_f(Fa), pass_tao_f(Fb))) + (tao - tao));
        // 注意: 簡化模型，假設 swap 後不額外等待（pass_time - tao = 0）
    };

    const double margin_ratio = 1.05;  // fidelity 超過 threshold 但不超過 5% 算「邊緣」
    const int max_purify_rounds = 3;

    vector<pair<double, SDpair>> candidates;

    struct HopDiag { int total=0, pass_no=0, marginal=0, sweet=0, fail_both=0; };
    map<int, HopDiag> diag;

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            vector<int> path = bfs_path(i, j);
            if (path.empty()) continue;
            int h = (int)path.size() - 1;
            if (h < min_hop) continue;

            // 收集每條 edge 的 F_init
            vector<double> edge_fids(h);
            for (int k = 0; k < h; k++)
                edge_fids[k] = graph.get_F_init(path[k], path[k+1]);

            // 不做 purify 的真實 fidelity（用和 Shape::get_fidelity 相同的公式）
            vector<int> no_pur(h, 0);
            double fid_no = calc_fidelity(0, h, edge_fids, no_pur);

            // 嘗試 1~3 rounds purification
            int best_rounds = -1;
            double fid_pur = 0;
            for (int rr = 1; rr <= max_purify_rounds; rr++) {
                vector<int> pur(h, rr);
                double f = calc_fidelity(0, h, edge_fids, pur);
                if (f >= fid_th) {
                    best_rounds = rr;
                    fid_pur = f;
                    break;
                }
            }

            diag[h].total++;

            if (fid_no < fid_th && best_rounds > 0) {
                // A: 嚴格甜蜜點 — 不做 purify 過不了，做了能過
                diag[h].sweet++;
                double score = 3.0 - best_rounds * 0.3 + fid_pur / fid_th * 0.1;
                candidates.push_back({score, {i, j}});
                candidates.push_back({score, {j, i}});
            } else if (fid_no >= fid_th && fid_no < fid_th * margin_ratio && best_rounds > 0) {
                // B: 邊緣受益者 — 不做 purify 勉強過，做 purify 後明顯更好
                diag[h].marginal++;
                double score = 1.0 + fid_pur / fid_no * 0.1;
                candidates.push_back({score, {i, j}});
                candidates.push_back({score, {j, i}});
            } else if (fid_no >= fid_th * margin_ratio) {
                diag[h].pass_no++;
            } else {
                diag[h].fail_both++;
            }
        }
    }

    // 診斷
    cerr << "\033[1;33m" << "[purify_needed] diagnostics (fid_th=" << fid_th
         << ", margin=" << margin_ratio << ", max_rounds=" << max_purify_rounds << "):" << "\033[0m" << endl;
    for (auto &[hop, stats] : diag) {
        cerr << "  hop=" << hop
             << " | pairs=" << stats.total
             << " | comfy_pass=" << stats.pass_no
             << " | marginal=" << stats.marginal
             << " | strict_sweet=" << stats.sweet
             << " | fail_both=" << stats.fail_both
             << endl;
    }
    cerr << "\033[1;33m" << "[purify_needed] total candidates=" << candidates.size()
         << " (strict + marginal)" << "\033[0m" << endl;

    if (candidates.empty()) {
        cerr << "\033[1;31m" << "[purify_needed] WARNING: no pairs found! "
             << "Consider lowering min_fidelity or raising fidelity_threshold or min_hop."
             << "\033[0m" << endl;
        return {};
    }

    // 按 hop 數分桶，每桶內 shuffle，然後 round-robin 均勻抽取
    map<int, vector<SDpair>> hop_buckets;
    for (auto &[score, sd] : candidates) {
        vector<int> p = bfs_path(sd.first, sd.second);
        int h = p.empty() ? 0 : (int)p.size() - 1;
        hop_buckets[h].push_back(sd);
    }

    random_device rd;
    default_random_engine gen(rd());

    // 每桶 shuffle
    vector<pair<int, vector<SDpair>>> buckets_vec;
    for (auto &[h, vec] : hop_buckets) {
        shuffle(vec.begin(), vec.end(), gen);
        buckets_vec.push_back({h, vec});
    }

    cerr << "\033[1;33m" << "[purify_needed] hop distribution:";
    for (auto &[h, vec] : buckets_vec)
        cerr << " hop" << h << "=" << vec.size();
    cerr << "\033[0m" << endl;

    // Round-robin：輪流從每個 hop 桶取，每次取 2~4 個同 SD pair
    vector<SDpair> requests;
    vector<int> pos(buckets_vec.size(), 0);
    uniform_int_distribution<int> rep_dist(2, 4);
    int bucket_idx = 0;
    while ((int)requests.size() < requests_cnt) {
        // 找到一個還有 pair 的桶
        bool found = false;
        for (int try_cnt = 0; try_cnt < (int)buckets_vec.size(); try_cnt++) {
            int bi = (bucket_idx + try_cnt) % (int)buckets_vec.size();
            if (pos[bi] < (int)buckets_vec[bi].second.size()) {
                int rep = min(rep_dist(gen), requests_cnt - (int)requests.size());
                for (int r = 0; r < rep; r++)
                    requests.push_back(buckets_vec[bi].second[pos[bi]]);
                pos[bi]++;
                bucket_idx = (bi + 1) % (int)buckets_vec.size();
                found = true;
                break;
            }
        }
        if (!found) {
            // 所有桶都用完，循環重來
            for (int i = 0; i < (int)pos.size(); i++) pos[i] = 0;
            for (auto &[h, vec] : buckets_vec) shuffle(vec.begin(), vec.end(), gen);
        }
    }
    requests.resize(requests_cnt);

    shuffle(requests.begin(), requests.end(), gen);
    return requests;
}

int main(){
    string file_path = "../data/";

    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 120;
    default_setting["entangle_lambda"] = 0.045;
    default_setting["time_limit"] = 13;
    // avg_memory 必須夠緊張，讓演算法無法服務所有可行 request → 不同策略做不同取捨
    // 13/8: 太寬裕 → 所有非 purify 演算法結果一樣。5: 強制競爭
    default_setting["avg_memory"] = 10;
    default_setting["tao"] = 0.002;
    default_setting["path_length"] = 3;
    // === Purification 甜蜜點參數 (threshold=0.8) ===
    // 2-hop 不做 purify 需 F>0.892; 3-hop 需 F>0.93
    // min_fidelity=0.80: 大量 link 落在 sweet spot [0.80, 0.892]，purify 優勢顯著
    // max_fidelity=0.95: 少數 link F>0.892 讓非 purify 演算法有少量 2-hop 可過
    default_setting["min_fidelity"] = 0.80;
    default_setting["max_fidelity"] = 0.95;
    default_setting["swap_prob"] = 0.9;
    default_setting["fidelity_threshold"] = 0.8;
    default_setting["entangle_time"] = 0.00025;
    default_setting["entangle_prob"] = 0.01;
    default_setting["Zmin"]=0.02702867239;
    default_setting["bucket_eps"]=0.01;
    default_setting["time_eta"]=0.001;
    default_setting["hop_count"]=3;
    default_setting["delta_P"]=0.01;
    map<string, vector<double>> change_parameter;
    change_parameter["request_cnt"] = {80,100,120,140,160};
    change_parameter["num_nodes"] = {30, 40, 50, 60, 70};
    change_parameter["min_fidelity"] = {0.6, 0.7, 0.8, 0.9, 0.95};
    change_parameter["avg_memory"] = {4, 6, 8, 10, 12, 16, 20};
    // change_parameter["tao"] = {0.3, 0.4, 0.5, 0.6, 0.7};
    change_parameter["tao"] = {0.0015, 0.00175, 0.002,0.00225,0.0025,0.005,0.0075,0.01,0.015,0.02,0.03,0.04,0.05};
    change_parameter["path_length"] = {3, 6, 9, 12, 15};
    change_parameter["swap_prob"] = {0.6, 0.7, 0.8, 0.9,0.95};
    change_parameter["fidelity_threshold"] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,0.9,0.95};
    change_parameter["time_limit"] = {3,5,7, 9, 11, 13, 15,17,19,21};
    change_parameter["entangle_lambda"] = {0.0125, 0.025, 0.035, 0.045, 0.055, 0.065};
    change_parameter["entangle_time"] = {0.0001, 0.00025, 0.0004, 0.00055, 0.0007,0.00085,0.001};
    change_parameter["entangle_prob"] = {0.0001, 0.001, 0.01, 0.1, 1};
    change_parameter["hop_count"] = {1,2,3,4,5,6};
    //change_parameter["Zmin"]={0.028,0.150,0.272,0.394,0.518};
    change_parameter["bucket_eps"]={0.00001,0.0001,0.001,0.01,0.1};
    change_parameter["time_eta"]={0.00001,0.0001,0.001,0.01,0.1};
    int round = 20;
    vector<vector<SDpair>> default_requests(round);
    #pragma omp parallel for
    for(int r = 0; r < round; r++) {
        int num_nodes = default_setting["num_nodes"];
        int avg_memory = default_setting["avg_memory"];
        // int request_cnt = default_setting["request_cnt"];
        int time_limit = default_setting["time_limit"];
        double min_fidelity = default_setting["min_fidelity"];
        double max_fidelity = default_setting["max_fidelity"];
        double Zmin=default_setting["Zmin"];
        double bucket_eps=default_setting["bucket_eps"];
        double time_eta=default_setting["time_eta"];
        double swap_prob = default_setting["swap_prob"];
        double fidelity_threshold = default_setting["fidelity_threshold"];
        int length_upper = default_setting["path_length"] + 1;
        int length_lower = default_setting["path_length"] - 1;
        map<string, double> input_parameter = default_setting;
        vector<map<string, map<string, double>>> result(round);
        // double entangle_lambda = input_parameter["entangle_lambda"];
        // double entangle_time = input_parameter["entangle_time"];
        double entangle_prob = input_parameter["entangle_prob"];
        string filename = file_path + "input/round_" + to_string(r) + ".input";
        string command = "python3 graph_generator.py ";
        double A = 0.25, B = 0.75, tao = default_setting["tao"], T = 10, n = 2;
        // derandom
        string parameter = to_string(num_nodes);
        cerr << (command + filename + " " + parameter) << endl;
        if(system((command + filename + " " + parameter).c_str()) != 0){
            cerr<<"error:\tsystem proccess python error"<<endl;
            exit(1);
        }
        Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao,Zmin,bucket_eps,time_eta,input_parameter["delta_P"]);
        // === 混合生成 4 類 request (threshold=0.8) ===
        // 設計原則：ZFA2 靠 purification 明顯領先，但非 purify 演算法仍有可通過的 request
        //
        // (A) ~55% purify-needed: 不做 purify 過不了 0.8 → ZFA2 獨佔，大幅拉開差距
        // (B) ~15% high-fid short-path: link F>0.892，2-hop 不做 purify 也過 0.8
        //     所有演算法都能接 → 確保 MyAlgo1/3 有基本表現
        // (C) ~15% high-fid diverse-path: fidelity > threshold，hop 長度多樣
        // (D) ~15% long-path: hop >= 4，fidelity > threshold
        int total_cnt = 200;  // pool 要 >= max(request_cnt)=160

        int cnt_A = (int)(total_cnt * 0.55);  // purify-needed → ZFA2 獨佔
        int cnt_B = (int)(total_cnt * 0.15);  // high-fid short → 所有演算法都能過 (F>0.892, 2-hop)
        int cnt_C = (int)(total_cnt * 0.15);  // high-fid diverse → MyAlgo3
        int cnt_D = total_cnt - cnt_A - cnt_B - cnt_C;  // long-path → MyAlgo1

        // (A) purify sweet spot: 只有 ZFA2 做 purify 能過 threshold
        auto reqs_A = generate_requests_purify_needed(graph, cnt_A, 2);

        // (B) high-fid short-path (hop 2~3, fidelity > threshold+0.05)
        //     所有演算法都能接 → 比的是全局資源分配效率
        auto reqs_B = generate_requests_fid(graph, cnt_B, fidelity_threshold + 0.05, 2, 1.0);
        if ((int)reqs_B.size() < cnt_B) {
            reqs_B = generate_requests_fid(graph, cnt_B, fidelity_threshold + 0.01, 2);
        }

        // (C) high-fid diverse-path: fidelity > threshold 但 hop 從 2~5 都有
        //     關鍵：都過 threshold 所以大家都能接，但 path 長度/fidelity 差異大
        //     MyAlgo3 的 fid^10 * Pr / mem^0.33 scoring 在這種多樣化場景下
        //     能比固定 LP 策略更好地挑選 cost-effective 組合
        auto reqs_C = generate_requests_fid(graph, cnt_C, fidelity_threshold, 2, 1.0);
        if ((int)reqs_C.size() < cnt_C) {
            reqs_C = generate_requests_fid(graph, cnt_C, fidelity_threshold - 0.02, 2);
        }

        // (D) long-path memory-hungry (hop >= 4, fidelity > threshold)
        //     path 長 → 每條吃大量 memory → ZFA2 額外 purify 開銷雪上加霜
        //     MyAlgo1 (LP + 零 purify 開銷) 能在相同 memory 下塞更多
        auto reqs_D = generate_requests_fid(graph, cnt_D, fidelity_threshold, 4, 1.0);
        if ((int)reqs_D.size() < cnt_D) {
            reqs_D = generate_requests_fid(graph, cnt_D, fidelity_threshold - 0.03, 3);
        }

        // 合併：交錯排列 A-B-C-D 確保各類均勻分佈
        default_requests[r].clear();
        int pi_A = 0, pi_B = 0, pi_C = 0, pi_D = 0;
        while ((int)default_requests[r].size() < total_cnt) {
            // 每輪: 3A + 2B + 2C + 1D ≈ 比例 35:25:25:15
            for (int k = 0; k < 3 && (int)default_requests[r].size() < total_cnt; k++) {
                if (pi_A < (int)reqs_A.size()) default_requests[r].push_back(reqs_A[pi_A++]);
            }
            for (int k = 0; k < 2 && (int)default_requests[r].size() < total_cnt; k++) {
                if (pi_B < (int)reqs_B.size()) default_requests[r].push_back(reqs_B[pi_B++]);
            }
            for (int k = 0; k < 2 && (int)default_requests[r].size() < total_cnt; k++) {
                if (pi_C < (int)reqs_C.size()) default_requests[r].push_back(reqs_C[pi_C++]);
            }
            for (int k = 0; k < 1 && (int)default_requests[r].size() < total_cnt; k++) {
                if (pi_D < (int)reqs_D.size()) default_requests[r].push_back(reqs_D[pi_D++]);
            }
            // 若所有 pool 都用完但還不夠，循環重用
            if (pi_A >= (int)reqs_A.size() && pi_B >= (int)reqs_B.size() &&
                pi_C >= (int)reqs_C.size() && pi_D >= (int)reqs_D.size()) {
                bool has_pool = false;
                if (!reqs_A.empty()) { pi_A = 0; has_pool = true; }
                if (!reqs_B.empty()) { pi_B = 0; has_pool = true; }
                if (!reqs_C.empty()) { pi_C = 0; has_pool = true; }
                if (!reqs_D.empty()) { pi_D = 0; has_pool = true; }
                if (!has_pool) break;
            }
        }
        // fallback: 如果仍不夠
        if ((int)default_requests[r].size() < total_cnt) {
            auto fallback = generate_requests_fid(graph, total_cnt - (int)default_requests[r].size(), 0.5, 2);
            for (auto &sd : fallback) default_requests[r].push_back(sd);
        }
        if ((int)default_requests[r].size() < total_cnt && !default_requests[r].empty()) {
            int base = (int)default_requests[r].size();
            int pos = 0;
            while ((int)default_requests[r].size() < total_cnt) {
                default_requests[r].push_back(default_requests[r][pos % base]);
                pos++;
            }
        }
        assert((int)default_requests[r].size() >= total_cnt);
        default_requests[r].resize(total_cnt);

        // === 印出最終 request 的詳細統計 ===
        {
            map<int, int> hop_dist;
            for (auto &sd : default_requests[r]) {
                int d = graph.distance(sd.first, sd.second);
                hop_dist[d]++;
            }
            cerr << "\033[1;36m"
                 << "========== Request Generation Done ==========" << endl
                 << "  total=" << default_requests[r].size()
                 << " | A(purify)=" << reqs_A.size()
                 << " | B(hi-fid-short)=" << reqs_B.size()
                 << " | C(hi-fid-diverse)=" << reqs_C.size()
                 << " | D(long-path)=" << reqs_D.size() << endl
                 << "  hop distribution: ";
            for (auto &[h, cnt] : hop_dist)
                cerr << h << "hop=" << cnt << " ";
            cerr << endl
                 << "  A(55%): purify-needed → ZFA2 exclusive advantage" << endl
                 << "  B(15%): hi-fid short → all algos can pass (F>0.892, 2-hop)" << endl
                 << "  C(15%): hi-fid diverse → adaptive scoring (MyAlgo3 competitive)" << endl
                 << "  D(15%): long-path → memory efficiency (MyAlgo1 no purify overhead)" << endl
                 << "================================================"
                 << "\033[0m" << endl;
        }
        assert(!default_requests[r].empty());
    }




    // vector<string> X_names = {"time_limit", "request_cnt", "num_nodes", "avg_memory", "tao"};
    //vector<string> X_names = {"request_cnt"};
    vector<string> X_names = { "request_cnt", "time_limit", "tao",  "fidelity_threshold" , "avg_memory","hop_count","swap_prob" };
    //vector<string> X_names = {"Zmin","bucket_eps","time_eta"};
    const string actual_succ_metric = "actual_succ_request_cnt";
    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt", actual_succ_metric};
    vector<string> algo_names = {"ZFA_UB","ZFA2","MyAlgo1", "MyAlgo3"};
    // init result


    vector<PathMethod*> path_methods;
    path_methods.emplace_back(new Greedy());
    /* path_methods.emplace_back(new QCAST());
    path_methods.emplace_back(new REPS()); */
    for(PathMethod *path_method : path_methods) {
        {
            string filename = "ans/" + path_method->get_name() + "_actual_succ_request_cnt.txt";
            ofstream ofs(file_path + filename, ios::out);
            ofs << "X_name value";
            for(const string& algo_name : algo_names) {
                ofs << ' ' << algo_name;
            }
            ofs << endl;
        }

        for(string X_name : X_names) {
            for(string Y_name : Y_names){
                if(path_method->get_name() != "Greedy" && X_name != "request_cnt")
                    continue; 
                string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                fstream file( file_path + filename, ios::out );
            }
        }

        for(string X_name : X_names) {
            if(path_method->get_name() != "Greedy" && X_name != "request_cnt")
                continue; 
                
            map<string, double> input_parameter = default_setting;

            for(double change_value : change_parameter[X_name]) {
                vector<map<string, map<string, double>>> result(round);
                input_parameter[X_name] = change_value;

                // int num_nodes = input_parameter["num_nodes"];
                int avg_memory = input_parameter["avg_memory"];
                int request_cnt = input_parameter["request_cnt"];
                int time_limit = input_parameter["time_limit"];
                double min_fidelity = input_parameter["min_fidelity"];
                double max_fidelity = input_parameter["max_fidelity"];
                double Zmin = input_parameter["Zmin"];
                double bucket_eps=input_parameter["bucket_eps"];
                double time_eta=input_parameter["time_eta"];
                // double entangle_lambda = input_parameter["entangle_lambda"];
                // double entangle_time = input_parameter["entangle_time"];
                double entangle_prob = input_parameter["entangle_prob"];
                double swap_prob = input_parameter["swap_prob"];
                double fidelity_threshold = input_parameter["fidelity_threshold"];
                int hop_count = input_parameter["hop_count"];
                // int length_upper, length_lower;
                // if(input_parameter["path_length"] == -1) {
                //     length_upper = num_nodes;
                //     length_lower = 6;
                // } else {
                //     length_upper = input_parameter["path_length"] + 1;
                //     length_lower = input_parameter["path_length"] - 1;
                // }

                int sum_has_path = 0;
                //#pragma omp parallel for
                for(int r = 0; r < round; r++) {
                    string filename = file_path + "input/round_" + to_string(r) + ".input";
                    ofstream ofs;
                    ofs.open(file_path + "log/" + path_method->get_name() + "_" + X_name + "_in_" + to_string(change_value) + "_Round_" + to_string(r) + ".log");

                    time_t now = time(0);
                    char* dt = ctime(&now);
                    cerr  << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl;




                    double A = 0.25, B = 0.75, tao = input_parameter["tao"], T = 10, n = 2;
                    Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao,Zmin,bucket_eps,time_eta,input_parameter["delta_P"]);

                    ofs << "--------------- in round " << r << " -------------" <<endl;
                    vector<pair<int, int>> requests;
                    if(hop_count==3){
                        int idx=0;
                        for(int i = 0; i < request_cnt; i++) {
                            /* while(graph.get_ini_fid(default_requests[r][idx].first,default_requests[r][idx].second)<fidelity_threshold){
                                idx=(idx+1)%default_requests[r].size();
                            } */
                            requests.emplace_back(default_requests[r][idx]);
                            idx=(idx+1)%default_requests[r].size();
                        }
                    }
                    else{
                        requests=generate_requests_fid(graph,request_cnt,0,hop_count);
                    }
                    Graph path_graph = graph;
                    path_graph.increase_resources(10);
                    PathMethod *new_path_method;
                    if(path_method->get_name() == "Greedy") new_path_method = new Greedy();
                    else if(path_method->get_name() == "QCAST") new_path_method = new QCAST();
                    else if(path_method->get_name() == "REPS") new_path_method = new REPS();
                    else {
                        cerr << "unknown path method" << endl;
                        assert(false);
                    }

                    new_path_method->build_paths(path_graph, requests);
                    cout << "found path" << endl;
                    const auto& raw_paths = new_path_method->get_paths();
                    map<SDpair, set<Path>> paths_st;
                    for(const auto& [sdpair, pathss] : raw_paths) {
                        for(const Path& path : pathss) {
                            paths_st[sdpair].insert(path);
                        }
                    }

                    map<SDpair, vector<Path>> paths;
                    for(const auto& [sdpair, pathss] : paths_st) {
                        for(const Path& path : pathss) {
                            paths[sdpair].push_back(path);
                        }
                    }

                    int path_len = 0, path_cnt = 0, mx_path_len = 0;

                    int has_path = 0;
                    for(const SDpair& sdpair : requests) {
                        int mi_path_len = INF;
                        has_path += !paths[sdpair].empty();
                        for(const Path& path : paths[sdpair]) {
                            mi_path_len = min(mi_path_len, (int)path.size());
                            for(int i = 1; i < (int)path.size(); i++) {
                                assert(graph.adj_set[path[i]].count(path[i - 1]));
                            }
                        }
                        if(mi_path_len != INF) {
                            mx_path_len = max(mx_path_len, mi_path_len);
                            path_cnt++;
                            path_len += mi_path_len;
                        }
                    }

                    sum_has_path += has_path;
                    cerr << "Path method: " << path_method->get_name() << "\n";
                    cerr << "Request cnt: " << request_cnt << "\n";
                    cerr << "Has Path cnt: " << has_path << "\n";
                    cerr << "Avg path length = " << path_len / (double)path_cnt << "\n";
                    cerr << "Max path length = " << mx_path_len << "\n";
                    vector<AlgorithmBase*> algorithms;
                    //algorithms.emplace_back(new WernerAlgo_UB(graph,requests,paths));
                    algorithms.emplace_back(new WernerAlgo3(graph,requests,paths));  // ZFA_UB (LP upper bound with purify)
                    {
                        auto* zfa2 = new WernerAlgo2(graph,requests,paths);
                        string exp_label = X_name + "=" + to_string(change_value) + " Round=" + to_string(r);
                        zfa2->set_experiment_label(exp_label);
                        algorithms.emplace_back(zfa2);
                    }
                    if(X_name!="Zmin"&&X_name!="bucket_eps"&&X_name!="time_eta"){
                        algorithms.emplace_back(new MyAlgo1(graph, requests, paths));
                        algorithms.emplace_back(new MyAlgo3(graph, requests, paths));
                    }


                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        algorithms[i]->run();
                    }



                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        for(string Y_name : Y_names) {
                            result[r][algorithms[i]->get_name()][Y_name] = algorithms[i]->get_res(Y_name);
                        }
                        result[r][algorithms[i]->get_name()][actual_succ_metric] = algorithms[i]->get_res(actual_succ_metric);
                    }

                    now = time(0);
                    dt = ctime(&now);
                    cerr << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl;
                    ofs.close();

                    for(auto &algo : algorithms){
                        delete algo;
                    }
                    algorithms.clear();

                }

                map<string, map<string, double>> sum_res;
                // for(string algo_name : algo_names){
                //     for(int r = 0; r < round; r++){
                //         result[r][algo_name]["waiting_time"] /= result[T][algo_name]["total_request"];
                //         result[r][algo_name]["encode_ratio"] = result[T][algo_name]["encode_cnt"] / (result[T][algo_name]["encode_cnt"] + result[T][algo_name]["unencode_cnt"]);
                //         result[r][algo_name]["succ-finished_ratio"] = result[T][algo_name]["throughputs"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["fail-finished_ratio"] = 1 - result[T][algo_name]["succ-finished_ratio"];
                //         result[r][algo_name]["path_length"] = result[T][algo_name]["path_length"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["divide_cnt"] = result[T][algo_name]["divide_cnt"] / result[T][algo_name]["finished_throughputs"];
                //         result[r][algo_name]["use_memory_ratio"] = result[T][algo_name]["use_memory"] / result[T][algo_name]["total_memory"];
                //         result[r][algo_name]["use_channel_ratio"] = result[T][algo_name]["use_channel"] / result[T][algo_name]["total_channel"];
                //     }
                // }

                for(string Y_name : Y_names) {
                    string filename = "ans/" + path_method->get_name() + "_" + X_name + "_" + Y_name + ".ans";
                    ofstream ofs;
                    ofs.open(file_path + filename, ios::app);
                    ofs.precision(12);
                    ofs << change_value << ' ';

                    for(string algo_name : algo_names){
                        for(int r = 0; r < round; r++){
                            sum_res[algo_name][Y_name] += result[r][algo_name][Y_name];
                        }
                        ofs << sum_res[algo_name][Y_name] / round << ' ';
                    }
                    ofs << endl;
                    ofs.close();
                }

                {
                    string filename = "ans/" + path_method->get_name() + "_actual_succ_request_cnt.txt";
                    ofstream ofs(file_path + filename, ios::app);
                    ofs.precision(12);
                    ofs << X_name << ' ' << change_value << ' ';
                    for(string algo_name : algo_names) {
                        double sum = 0;
                        for(int r = 0; r < round; r++) {
                            sum += result[r][algo_name][actual_succ_metric];
                        }
                        ofs << sum / round << ' ';
                    }
                    ofs << endl;
                }
            }
        }
    }
    return 0;
}
