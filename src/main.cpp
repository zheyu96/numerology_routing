
#include "./config.h"
#include <sys/resource.h>
#include <new>
#include <stdexcept>
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/MyAlgo1/MyAlgo1.h"
#include "Algorithm/MyAlgo2/MyAlgo2.h"
#include "Algorithm/MyAlgo3/MyAlgo3.h"
#include "Algorithm/MyAlgo4/MyAlgo4.h"
#include "Algorithm/MyAlgo5/MyAlgo5.h"
#include "Algorithm/MyAlgo6/MyAlgo6.h"
#include "Algorithm/GreedySP/GreedySP.h"
// PathMethod 系列必須放在 Werner 系列之前：Werner*.h 內有
// `#define double long double`，會污染其後 include 的所有 header
#include "Network/PathMethod/PathMethodBase/PathMethod.h"
#include "Network/PathMethod/Greedy/Greedy.h"
#include "Network/PathMethod/QCAST/QCAST.h"
#include "Network/PathMethod/REPS/REPS.h"
#include "Algorithm/WernerAlgo/WernerAlgo.h"
#include "Algorithm/WernerAlgo2/WernerAlgo2.h"
#include "Algorithm/WernerAlgo_routing/WernerAlgo_routing.h"
#include "Algorithm/WernerAlgo3/WernerAlgo3.h"
#include "Algorithm/WernerAlgo_UB/WernerAlgo_UB.h"

using namespace std;

// [DEBUG] 印出當前 RSS / 系統記憶體 — 方便定位 std::bad_alloc 發生點
static void DBG_mem(const char* tag) {
    struct rusage ru;
    if(getrusage(RUSAGE_SELF, &ru) == 0) {
        // ru_maxrss 在 Linux 是 KB
        cerr << "[MEM] " << tag << " RSS=" << ru.ru_maxrss << " KB" << endl;
    } else {
        cerr << "[MEM] " << tag << " (getrusage failed)" << endl;
    }
    cerr.flush();
}

#define DBG_HERE(tag) do { cerr << "[CKPT] " << tag << endl; cerr.flush(); } while(0)

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

// A tau sweep must use the same workload at every x-axis point.  Build that
// workload once, but classify every source/destination pair against the whole
// sweep instead of filtering at only the default tau.
enum class RequestTier {
    ROBUST,             // feasible at the largest tau
    MID_ONLY,           // feasible at the median tau, not at the largest tau
    LOW_ONLY            // feasible at the smallest tau only
};

struct RequestCandidate {
    SDpair sd;
    int hop = 0;
    RequestTier tier = RequestTier::LOW_ONLY;
};

static Path shortest_path(Graph& graph, int src, int dst) {
    const int n = graph.get_num_nodes();
    vector<int> parent(n, -1);
    queue<int> q;
    parent[src] = src;
    q.push(src);

    while(!q.empty()) {
        int u = q.front();
        q.pop();
        if(u == dst) break;
        for(int v : graph.adj_list[u]) {
            if(parent[v] != -1) continue;
            parent[v] = u;
            q.push(v);
        }
    }
    if(parent[dst] == -1) return {};

    Path path;
    for(int v = dst; v != src; v = parent[v]) path.push_back(v);
    path.push_back(src);
    reverse(path.begin(), path.end());
    return path;
}

// Fidelity estimator matching Shape::recursion_get_fidelity for the balanced
// schedule used by the request filters.  With purify_rounds == 0 every leaf
// starts from the raw elementary-link fidelity; purify_rounds > 0 pumps every
// leaf that many times (paper Eq.(4), same model as the legacy purify-aware
// generator) before decoherence.  Purification's extra time slots are not
// modelled here, so keep purify_rounds small for honest estimates.
static double estimate_balanced_fidelity(Graph& graph, const Path& path,
                                         double tau, int purify_rounds = 0) {
    const double A = graph.get_A();
    const double B = graph.get_B();
    const double T = graph.get_T();
    const double exponent = graph.get_n();

    auto t2F = [&](double t) -> double {
        if(t >= 1e5) return 0.0;
        return A + B * exp(-pow(t / T, exponent));
    };
    auto F2t = [&](double F) -> double {
        if(F <= A + EPS) return 1e9;
        return T * pow(-log((F - A) / B), 1.0 / exponent);
    };
    auto pass_tau = [&](double F) -> double {
        return t2F(F2t(F) + tau);
    };
    auto swap_fidelity = [&](double left, double right) -> double {
        if(left <= A + EPS || right <= A + EPS) return 0.0;
        return left * right + (1.0 / 3.0) * (1.0 - left) * (1.0 - right);
    };

    auto pump = [&](double F) -> double {
        double Fbase = F, Fcur = F;
        for(int r = 0; r < purify_rounds; r++) {
            double num = Fcur * Fbase
                       + (1.0 / 9.0) * (1.0 - Fcur) * (1.0 - Fbase);
            double den = Fcur * Fbase
                       + (1.0 / 3.0) * Fcur * (1.0 - Fbase)
                       + (1.0 / 3.0) * (1.0 - Fcur) * Fbase
                       + (5.0 / 9.0) * (1.0 - Fcur) * (1.0 - Fbase);
            Fcur = num / den;
        }
        return Fcur;
    };

    function<double(int, int)> solve = [&](int left, int right) -> double {
        if(left + 1 == right) {
            double fidelity = pump(graph.get_F_init(path[left], path[right]));
            return pass_tau(fidelity);
        }
        int mid = (left + right) / 2;
        double left_fidelity = solve(left, mid);
        double right_fidelity = solve(mid, right);
        return swap_fidelity(pass_tau(left_fidelity), pass_tau(right_fidelity));
    };

    return solve(0, (int)path.size() - 1);
}

static vector<SDpair> generate_tau_stratified_requests(
        Graph& graph, int request_count, vector<double> tau_values,
        unsigned int seed) {
    if(request_count <= 0 || tau_values.empty()) return {};
    sort(tau_values.begin(), tau_values.end());
    tau_values.erase(unique(tau_values.begin(), tau_values.end()), tau_values.end());

    const double threshold = graph.get_fidelity_threshold();
    const double tau_low = tau_values.front();
    const double tau_mid = tau_values[tau_values.size() / 2];
    const double tau_high = tau_values.back();
    const int n = graph.get_num_nodes();

    map<RequestTier, vector<RequestCandidate>> buckets;
    map<int, int> all_hops;
    int rejected_at_low_tau = 0;

    // Use unordered node pairs for classification, then add both directions.
    // Direction is symmetric physically but keeping both prevents orientation
    // bias in algorithms whose tie-breaking depends on node order.
    for(int src = 0; src < n; ++src) {
        for(int dst = src + 1; dst < n; ++dst) {
            Path path = shortest_path(graph, src, dst);
            if(path.size() < 2) continue;

            int hop = (int)path.size() - 1;
            double high_fidelity = estimate_balanced_fidelity(graph, path, tau_high);
            double mid_fidelity = estimate_balanced_fidelity(graph, path, tau_mid);
            double low_fidelity = estimate_balanced_fidelity(graph, path, tau_low);

            RequestTier tier;
            if(high_fidelity + EPS >= threshold)
                tier = RequestTier::ROBUST;
            else if(mid_fidelity + EPS >= threshold)
                tier = RequestTier::MID_ONLY;
            else if(low_fidelity + EPS >= threshold)
                tier = RequestTier::LOW_ONLY;
            else {
                rejected_at_low_tau++;
                continue;
            }

            buckets[tier].push_back({{src, dst}, hop, tier});
            buckets[tier].push_back({{dst, src}, hop, tier});
            all_hops[hop] += 2;
        }
    }

    mt19937 rng(seed);
    for(auto& [tier, candidates] : buckets)
        shuffle(candidates.begin(), candidates.end(), rng);

    auto tier_name = [](RequestTier tier) {
        switch(tier) {
            case RequestTier::ROBUST: return "robust";
            case RequestTier::MID_ONLY: return "mid-only";
            case RequestTier::LOW_ONLY: return "low-only";
        }
        return "unknown";
    };

    cerr << "[request-filter] tau range=" << tau_low << ".." << tau_high
         << " mid=" << tau_mid << " threshold=" << threshold
         << " purification=disabled" << endl;
    for(RequestTier tier : {RequestTier::ROBUST, RequestTier::MID_ONLY,
                            RequestTier::LOW_ONLY}) {
        cerr << "  " << tier_name(tier) << " candidates=" << buckets[tier].size() << endl;
    }
    cerr << "  rejected-even-at-low-tau pairs=" << rejected_at_low_tau * 2 << endl;

    if(buckets[RequestTier::ROBUST].empty()) {
        cerr << "[request-filter] WARNING: no request is feasible at maximum tau="
             << tau_high << "; a non-zero high-tau result cannot be guaranteed with "
             << "the current physical parameters." << endl;
    }

    // Half of the workload is feasible at tau_high without purification.  The
    // other half preserves medium/low-tau stress cases so the curve still
    // measures degradation rather than becoming an all-easy benchmark.
    map<RequestTier, int> target;
    target[RequestTier::ROBUST] = request_count * 50 / 100;
    target[RequestTier::MID_ONLY] = request_count * 30 / 100;
    target[RequestTier::LOW_ONLY] = request_count
        - target[RequestTier::ROBUST]
        - target[RequestTier::MID_ONLY];

    vector<SDpair> requests;
    map<RequestTier, size_t> cursor;
    map<RequestTier, int> selected_by_tier;

    auto append_from = [&](RequestTier preferred, int count,
                           initializer_list<RequestTier> fallback_order) {
        for(int k = 0; k < count; ++k) {
            RequestTier chosen = preferred;
            bool found = false;
            for(RequestTier candidate_tier : fallback_order) {
                if(!buckets[candidate_tier].empty()) {
                    chosen = candidate_tier;
                    found = true;
                    break;
                }
            }
            if(!found) return;

            auto& candidates = buckets[chosen];
            if(cursor[chosen] > 0 && cursor[chosen] % candidates.size() == 0)
                shuffle(candidates.begin(), candidates.end(), rng);
            const auto& candidate = candidates[cursor[chosen] % candidates.size()];
            cursor[chosen]++;
            requests.push_back(candidate.sd);
            selected_by_tier[chosen]++;
        }
    };

    // Interleave each group of 20 requests (10/6/4).  Consequently the common
    // request-count prefixes 80/100/120/... preserve the intended mix and the
    // 50% robust high-tau guarantee.
    vector<RequestTier> tier_schedule;
    while((int)tier_schedule.size() + 20 <= request_count) {
        vector<RequestTier> cycle;
        cycle.insert(cycle.end(), 10, RequestTier::ROBUST);
        cycle.insert(cycle.end(), 6, RequestTier::MID_ONLY);
        cycle.insert(cycle.end(), 4, RequestTier::LOW_ONLY);
        shuffle(cycle.begin(), cycle.end(), rng);
        tier_schedule.insert(tier_schedule.end(), cycle.begin(), cycle.end());
    }
    for(RequestTier tier : {RequestTier::ROBUST, RequestTier::MID_ONLY,
                            RequestTier::LOW_ONLY}) {
        int already = count(tier_schedule.begin(), tier_schedule.end(), tier);
        int missing = max(0, target[tier] - already);
        tier_schedule.insert(tier_schedule.end(), missing, tier);
    }
    shuffle(tier_schedule.begin() + (tier_schedule.size() / 20) * 20,
            tier_schedule.end(), rng);

    for(RequestTier tier : tier_schedule) {
        if(tier == RequestTier::ROBUST)
            append_from(tier, 1, {RequestTier::ROBUST, RequestTier::MID_ONLY,
                                  RequestTier::LOW_ONLY});
        else if(tier == RequestTier::MID_ONLY)
            append_from(tier, 1, {RequestTier::MID_ONLY, RequestTier::ROBUST,
                                  RequestTier::LOW_ONLY});
        else
            append_from(tier, 1, {RequestTier::LOW_ONLY, RequestTier::MID_ONLY,
                                  RequestTier::ROBUST});
    }

    if((int)requests.size() != request_count) {
        throw runtime_error("request filter found no pair feasible even at minimum tau");
    }
    map<int, int> selected_hops;
    map<double, int> predicted_feasible;
    for(const SDpair& sd : requests) {
        Path path = shortest_path(graph, sd.first, sd.second);
        selected_hops[(int)path.size() - 1]++;
        for(double tau : tau_values) {
            double fidelity = estimate_balanced_fidelity(graph, path, tau);
            if(fidelity + EPS >= threshold) predicted_feasible[tau]++;
        }
    }

    cerr << "[request-filter] selected=" << requests.size() << " tiers:";
    for(RequestTier tier : {RequestTier::ROBUST, RequestTier::MID_ONLY,
                            RequestTier::LOW_ONLY})
        cerr << " " << tier_name(tier) << "=" << selected_by_tier[tier];
    cerr << endl << "[request-filter] selected hop distribution:";
    for(auto [hop, count] : selected_hops) cerr << " " << hop << "hop=" << count;
    cerr << endl << "[request-filter] predicted feasible (before competition):";
    for(auto [tau, count] : predicted_feasible)
        cerr << " tau=" << tau << ":" << count << "/" << requests.size();
    cerr << endl;

    return requests;
}

// ===== Routing-favored workload（能力階梯過濾器）=====
// 依「演算法能力階梯」分類每個 SD pair，讓每一級只有能力更強的演算法能服務：
//   SP_OK       : 最短路徑(不 purify)在整個 tau sweep 都達標
//                 → 所有演算法都能服務（維持 baseline 基本盤，benchmark 才有鑑別度）
//   NEED_PURIFY : 最短路徑不 purify 在 tau_low 就過不了（fidelity 隨 tau 單調下降
//                 ⇒ 整個 sweep 都過不了），purify(≤MAX_PURIFY_ROUNDS) 後 tau_mid 可過
//                 → 只有會 purify 的 ZFA 系列（ZFA_UB / ZFA2 / ZFA_routing）能服務
//   NEED_DETOUR : 所有最短路徑連 purify 都救不回，但存在 +DETOUR_SLACK hop 的
//                 替代路徑不 purify 就在 tau_mid 達標
//                 → path set 只含最短路徑，所以只有自由選路的 ZFA_routing 能服務
//   NEED_BOTH   : 替代路徑也要 purify 才在 tau_mid 達標
//                 → 只有 ZFA_routing 能服務，且必須 purify
// 判定點：「過不了」用 tau_low（單調性 ⇒ 整個 sweep 都過不了）；「過得了」用
// tau_mid（⇒ tau ≤ tau_mid 的實驗點保證可服務；tau > tau_mid 各類別自然退化）。
enum class WorkloadClass { SP_OK, NEED_PURIFY, NEED_DETOUR, NEED_BOTH };

struct RoutingFavoredCandidate {
    SDpair sd;
    int sp_hop = 0;
    WorkloadClass cls = WorkloadClass::SP_OK;
    // 各值皆為「該類所有候選路徑的最大 fidelity」，供分類與診斷用
    double sp0_low = 0, sp0_mid = 0, sp0_high = 0;   // 最短路徑、不 purify
    double spP_low = 0, spP_mid = 0, spP_high = 0;   // 最短路徑、purify 後
    double alt0_mid = 0, alt0_high = 0;              // 替代路徑、不 purify
    double altP_mid = 0, altP_high = 0;              // 替代路徑、purify 後
};

static vector<SDpair> generate_routing_favored_requests(
        Graph& graph, int request_count, vector<double> tau_values,
        unsigned int seed) {
    if(request_count <= 0 || tau_values.empty()) return {};
    sort(tau_values.begin(), tau_values.end());
    tau_values.erase(unique(tau_values.begin(), tau_values.end()), tau_values.end());

    const double threshold = graph.get_fidelity_threshold();
    const double tau_low = tau_values.front();
    const double tau_mid = tau_values[tau_values.size() / 2];
    const double tau_high = tau_values.back();
    const int n = graph.get_num_nodes();

    const int MAX_SP_HOP = 4;         // 最短路徑 hop 上限（time_limit 內可完成）
    const int DETOUR_SLACK = 1;       // 替代路徑允許比最短路徑多的 hop 數
    const int MAX_PURIFY_ROUNDS = 5;  // 與 WPFA/SP 實際演算法一致
    const int TOPK = 12;              // 保留 high-tau 最佳路徑做 purify 評估
    const int MAX_PATHS_PER_PAIR = 20000;

    auto bfs_dist = [&](int root) {
        vector<int> dist(n, INF);
        queue<int> que;
        dist[root] = 0;
        que.push(root);
        while(!que.empty()) {
            int u = que.front(); que.pop();
            for(int v : graph.adj_list[u]) {
                if(dist[v] > dist[u] + 1) {
                    dist[v] = dist[u] + 1;
                    que.push(v);
                }
            }
        }
        return dist;
    };

    map<WorkloadClass, vector<RoutingFavoredCandidate>> buckets;
    int rejected_infeasible = 0;   // 連 ZFA_routing 都救不回
    int rejected_tau_marginal = 0; // 不 purify 在 tau_low 可過但 tau_high 不行（模糊帶）

    for(int dst = 0; dst < n; ++dst) {
        vector<int> dist_to_dst = bfs_dist(dst);
        for(int src = 0; src < dst; ++src) {
            int sp = dist_to_dst[src];
            if(sp < 1 || sp > MAX_SP_HOP) continue;
            const int max_hop = sp + DETOUR_SLACK;

            RoutingFavoredCandidate cand;
            cand.sd = {src, dst};
            cand.sp_hop = sp;
            vector<pair<double, Path>> top_sp, top_alt;  // (f_high, path)
            int path_cnt = 0;

            auto keep_topk = [&](vector<pair<double, Path>>& top,
                                 double f, const Path& p) {
                if((int)top.size() < TOPK) {
                    top.push_back({f, p});
                    return;
                }
                int worst = 0;
                for(int i = 1; i < (int)top.size(); i++)
                    if(top[i].first < top[worst].first) worst = i;
                if(f > top[worst].first) top[worst] = {f, p};
            };

            // DFS 枚舉 hop ≤ max_hop 的所有 simple path（dist_to_dst 剪枝）
            Path cur = {src};
            vector<bool> vis(n, false);
            vis[src] = true;
            function<void(int)> dfs = [&](int u) {
                if(path_cnt >= MAX_PATHS_PER_PAIR) return;
                if(u == dst) {
                    path_cnt++;
                    int h = (int)cur.size() - 1;
                    double f_mid = estimate_balanced_fidelity(graph, cur, tau_mid);
                    double f_high = estimate_balanced_fidelity(graph, cur, tau_high);
                    if(h == sp) {
                        cand.sp0_mid = max(cand.sp0_mid, f_mid);
                        cand.sp0_low = max(cand.sp0_low,
                            estimate_balanced_fidelity(graph, cur, tau_low));
                        cand.sp0_high = max(cand.sp0_high, f_high);
                        keep_topk(top_sp, f_high, cur);
                    } else {
                        cand.alt0_mid = max(cand.alt0_mid, f_mid);
                        cand.alt0_high = max(cand.alt0_high, f_high);
                        keep_topk(top_alt, f_high, cur);
                    }
                    return;
                }
                for(int v : graph.adj_list[u]) {
                    if(vis[v]) continue;
                    if((int)cur.size() + dist_to_dst[v] > max_hop) continue;
                    vis[v] = true;
                    cur.push_back(v);
                    dfs(v);
                    cur.pop_back();
                    vis[v] = false;
                }
            };
            dfs(src);
            if(top_sp.empty()) continue;

            // purify 評估只做在 top-K 路徑上（以 f_mid 排序的近似；pumping 對
            // 每條 link 單調，排名通常一致，偏差只影響 workload 組成不影響正確性）
            for(const auto& [f0, p] : top_sp) {
                for(int rr = 1; rr <= MAX_PURIFY_ROUNDS; rr++) {
                    cand.spP_mid = max(cand.spP_mid,
                        estimate_balanced_fidelity(graph, p, tau_mid, rr));
                    cand.spP_low = max(cand.spP_low,
                        estimate_balanced_fidelity(graph, p, tau_low, rr));
                    cand.spP_high = max(cand.spP_high,
                        estimate_balanced_fidelity(graph, p, tau_high, rr));
                }
            }
            for(const auto& [f0, p] : top_alt) {
                for(int rr = 1; rr <= MAX_PURIFY_ROUNDS; rr++) {
                    cand.altP_mid = max(cand.altP_mid,
                        estimate_balanced_fidelity(graph, p, tau_mid, rr));
                    cand.altP_high = max(cand.altP_high,
                        estimate_balanced_fidelity(graph, p, tau_high, rr));
                }
            }

            WorkloadClass cls;
            if(cand.sp0_high + EPS >= threshold) {
                cls = WorkloadClass::SP_OK;
            } else if(cand.sp0_low + EPS < threshold
                      && cand.spP_high + EPS >= threshold) {
                cls = WorkloadClass::NEED_PURIFY;
            } else if(cand.spP_low + EPS < threshold
                      && cand.alt0_high + EPS >= threshold) {
                cls = WorkloadClass::NEED_DETOUR;
            } else if(cand.spP_low + EPS < threshold
                      && cand.altP_high + EPS >= threshold) {
                cls = WorkloadClass::NEED_BOTH;
            } else if(cand.sp0_low + EPS >= threshold) {
                rejected_tau_marginal++;
                continue;
            } else {
                rejected_infeasible++;
                continue;
            }

            cand.cls = cls;
            buckets[cls].push_back(cand);
            RoutingFavoredCandidate rev = cand;
            swap(rev.sd.first, rev.sd.second);
            buckets[cls].push_back(rev);
        }
    }

    auto cls_name = [](WorkloadClass c) {
        switch(c) {
            case WorkloadClass::SP_OK:       return "sp-ok";
            case WorkloadClass::NEED_PURIFY: return "need-purify";
            case WorkloadClass::NEED_DETOUR: return "need-detour";
            case WorkloadClass::NEED_BOTH:   return "need-both";
        }
        return "unknown";
    };
    const vector<WorkloadClass> ALL_CLS = {
        WorkloadClass::SP_OK, WorkloadClass::NEED_PURIFY,
        WorkloadClass::NEED_DETOUR, WorkloadClass::NEED_BOTH};

    cerr << "[routing-favored] tau low/mid/high=" << tau_low << "/" << tau_mid
         << "/" << tau_high << " threshold=" << threshold
         << " detour_slack=" << DETOUR_SLACK
         << " max_purify_rounds=" << MAX_PURIFY_ROUNDS << endl;
    for(WorkloadClass c : ALL_CLS)
        cerr << "  " << cls_name(c) << " candidates=" << buckets[c].size() << endl;
    cerr << "  rejected: infeasible=" << rejected_infeasible * 2
         << " tau-marginal=" << rejected_tau_marginal * 2 << endl;
    if(buckets[WorkloadClass::NEED_DETOUR].empty()
       && buckets[WorkloadClass::NEED_BOTH].empty()) {
        cerr << "[routing-favored] WARNING: no detour-advantaged pair exists; "
             << "ZFA_routing has no exclusive request.  Widen the link-quality "
             << "spread (e.g. lower min_fidelity) to create them." << endl;
    }

    mt19937 rng(seed);
    for(auto& [cls, cands] : buckets)
        shuffle(cands.begin(), cands.end(), rng);

    // 每 20 個 request 中，14 個需要 purification；其中 6 個還需要 detour。
    // 另外保留 4 個 SP_OK 與 2 個 detour-only 作公平對照。
    vector<WorkloadClass> schedule;
    while((int)schedule.size() < request_count) {
        vector<WorkloadClass> cycle;
        cycle.insert(cycle.end(), 4, WorkloadClass::SP_OK);
        cycle.insert(cycle.end(), 8, WorkloadClass::NEED_PURIFY);
        cycle.insert(cycle.end(), 2, WorkloadClass::NEED_DETOUR);
        cycle.insert(cycle.end(), 6, WorkloadClass::NEED_BOTH);
        shuffle(cycle.begin(), cycle.end(), rng);
        schedule.insert(schedule.end(), cycle.begin(), cycle.end());
    }
    schedule.resize(request_count);

    vector<SDpair> requests;
    vector<RoutingFavoredCandidate> selected;
    map<WorkloadClass, size_t> cursor;
    map<WorkloadClass, int> selected_by_cls;

    auto take_one = [&](initializer_list<WorkloadClass> order) {
        for(WorkloadClass c : order) {
            auto& cands = buckets[c];
            if(cands.empty()) continue;
            if(cursor[c] > 0 && cursor[c] % cands.size() == 0)
                shuffle(cands.begin(), cands.end(), rng);
            const RoutingFavoredCandidate cand = cands[cursor[c] % cands.size()];
            cursor[c]++;
            requests.push_back(cand.sd);
            selected.push_back(cand);
            selected_by_cls[c]++;
            return true;
        }
        return false;
    };

    for(WorkloadClass c : schedule) {
        bool ok = false;
        // fallback 順序：先往「同樣偏好 routing」的類別找，最後才退到 sp-ok
        switch(c) {
            case WorkloadClass::SP_OK:
                ok = take_one({WorkloadClass::SP_OK, WorkloadClass::NEED_PURIFY,
                               WorkloadClass::NEED_DETOUR, WorkloadClass::NEED_BOTH});
                break;
            case WorkloadClass::NEED_PURIFY:
                ok = take_one({WorkloadClass::NEED_PURIFY, WorkloadClass::NEED_BOTH,
                               WorkloadClass::NEED_DETOUR, WorkloadClass::SP_OK});
                break;
            case WorkloadClass::NEED_DETOUR:
                ok = take_one({WorkloadClass::NEED_DETOUR, WorkloadClass::NEED_BOTH,
                               WorkloadClass::NEED_PURIFY, WorkloadClass::SP_OK});
                break;
            case WorkloadClass::NEED_BOTH:
                ok = take_one({WorkloadClass::NEED_BOTH, WorkloadClass::NEED_DETOUR,
                               WorkloadClass::NEED_PURIFY, WorkloadClass::SP_OK});
                break;
        }
        if(!ok)
            throw runtime_error("routing-favored filter: every bucket is empty");
    }

    map<int, int> hop_dist;
    for(const auto& c : selected) hop_dist[c.sp_hop]++;
    // 預測各能力等級在 tau_mid 可服務的 request 數（尚未考慮資源競爭）
    int base_low = 0, base_mid = 0, base_high = 0;
    int zfa_mid = 0, zfa_high = 0, routing_mid = 0, routing_high = 0;
    for(const auto& c : selected) {
        if(c.sp0_low + EPS >= threshold) base_low++;
        if(c.sp0_mid + EPS >= threshold) base_mid++;
        if(c.sp0_high + EPS >= threshold) base_high++;
        if(max(c.sp0_mid, c.spP_mid) + EPS >= threshold) zfa_mid++;
        if(max(c.sp0_high, c.spP_high) + EPS >= threshold) zfa_high++;
        if(max(max(c.sp0_mid, c.spP_mid),
               max(c.alt0_mid, c.altP_mid)) + EPS >= threshold) routing_mid++;
        if(max(max(c.sp0_high, c.spP_high),
               max(c.alt0_high, c.altP_high)) + EPS >= threshold) routing_high++;
    }

    cerr << "[routing-favored] selected=" << requests.size() << " mix:";
    for(WorkloadClass c : ALL_CLS)
        cerr << " " << cls_name(c) << "=" << selected_by_cls[c];
    cerr << endl << "[routing-favored] sp-hop distribution:";
    for(auto [hop, count] : hop_dist) cerr << " " << hop << "hop=" << count;
    cerr << endl << "[routing-favored] predicted servable (before competition):"
         << " @tau_mid SP-no-purify=" << base_mid
         << " ZFA(purify)=" << zfa_mid
         << " ZFA_routing=" << routing_mid << "/" << requests.size()
         << " | @tau_high ZFA(purify)=" << zfa_high
         << " ZFA_routing=" << routing_high
         << " | SP-no-purify @tau_low=" << base_low
         << " @tau_high=" << base_high << endl;

    return requests;
}
int main(){
    string file_path = "../data/";

    map<string, double> default_setting;
    // Small, reproducible default for quick tau validation.  Larger runs can
    // override these with WPFA_NUM_NODES / WPFA_REQUEST_COUNT.
    default_setting["num_nodes"] = 30;
    default_setting["request_cnt"] = 30;
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
    // Paper TRIM_epsilon is tunable. 0.01 makes the all-pairs routing
    // extension impractically large after removing the non-paper hard caps.
    default_setting["bucket_eps"]=0.05;
    default_setting["time_eta"]=0.001;
    default_setting["hop_count"]=3;
    default_setting["delta_P"]=0.01;
    auto apply_int_env = [&](const char* name, const char* key) {
        if(const char* value = getenv(name)) {
            try {
                default_setting[key] = max(1, stoi(value));
            } catch(const exception&) {
                cerr << "[config] ignoring invalid " << name << "='" << value << "'" << endl;
            }
        }
    };
    apply_int_env("WPFA_NUM_NODES", "num_nodes");
    apply_int_env("WPFA_REQUEST_COUNT", "request_cnt");
    map<string, vector<double>> change_parameter;
    change_parameter["request_cnt"] = {80,100,120,140,160};
    change_parameter["num_nodes"] = {30, 40, 50, 60, 70};
    change_parameter["min_fidelity"] = {0.6, 0.7, 0.8, 0.9, 0.95};
    change_parameter["avg_memory"] = {4, 6, 8, 10, 12, 16, 20};
    change_parameter["tao"] = {0.01, 0.02, 0.03, 0.04, 0.05};
    change_parameter["path_length"] = {3, 6, 9, 12, 15};
    change_parameter["swap_prob"] = {0.6, 0.7, 0.8, 0.9,0.95};
    change_parameter["fidelity_threshold"] = {0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,0.9,0.95};
    change_parameter["time_limit"] = {5,7, 9, 11, 13, 15,17,19};
    change_parameter["entangle_lambda"] = {0.0125, 0.025, 0.035, 0.045, 0.055, 0.065};
    change_parameter["entangle_time"] = {0.0001, 0.00025, 0.0004, 0.00055, 0.0007,0.00085,0.001};
    change_parameter["entangle_prob"] = {0.0001, 0.001, 0.01, 0.1, 1};
    change_parameter["hop_count"] = {1,2,3,4,5,6};
    //change_parameter["Zmin"]={0.028,0.150,0.272,0.394,0.518};
    change_parameter["bucket_eps"]={0.00001,0.0001,0.001,0.01,0.1};
    change_parameter["time_eta"]={0.00001,0.0001,0.001,0.01,0.1};
    int round = 5;
    if(const char* env_rounds = getenv("WPFA_ROUNDS")) {
        try {
            round = max(1, stoi(env_rounds));
        } catch(const exception&) {
            cerr << "[config] invalid WPFA_ROUNDS='" << env_rounds
                 << "'; using 1" << endl;
            round = 1;
        }
    }
    // Run the complete algorithm set by default, even with the small quick
    // parameters.  WPFA_COMPARE_ONLY=1 is an optional diagnostic shortcut.
    const bool compare_only = getenv("WPFA_COMPARE_ONLY") != nullptr;
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
#ifdef _WIN32
        string command = "python graph_generator.py ";
#else
        string command = "python3 graph_generator.py ";
#endif
        double A = 0.25, B = 0.75, tao = default_setting["tao"], T = 10, n = 2;
        // derandom
        string parameter = to_string(num_nodes);
        cerr << (command + filename + " " + parameter) << endl;
        if(system((command + filename + " " + parameter).c_str()) != 0){
            cerr<<"error:\tsystem proccess python error"<<endl;
            exit(1);
        }
        Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao,Zmin,bucket_eps,time_eta,input_parameter["delta_P"],input_parameter["entangle_lambda"],input_parameter["entangle_time"]);
#if 0  // Legacy purification-aware request generator (disabled)
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
#endif

        // Build one deterministic workload per topology and classify every
        // candidate against the complete tau sweep.  Every experiment point
        // later reuses this exact request sequence.
        // This is only an offline candidate pool.  It does not change the
        // online request_cnt passed to any algorithm.
        const int request_pool_size = 200;
        // true : 能力階梯 workload（sp-ok / need-purify / need-detour / need-both）
        // false: 原本的 tau 分層 workload（不 purify、只看最短路徑）
        const bool routing_favored_workload = true;
        default_requests[r] = routing_favored_workload
            ? generate_routing_favored_requests(
                  graph, request_pool_size, change_parameter["tao"],
                  0x5EEDu + static_cast<unsigned int>(r))
            : generate_tau_stratified_requests(
                  graph, request_pool_size, change_parameter["tao"],
                  0x5EEDu + static_cast<unsigned int>(r));

        map<int, int> hop_dist;
        for(const auto& sd : default_requests[r])
            hop_dist[graph.distance(sd.first, sd.second)]++;
        cerr << "========== Request Generation Done ==========\n"
             << "  total=" << default_requests[r].size()
             << " | generator="
             << (routing_favored_workload ? "routing-favored-capability-ladder"
                                          : "tau-stratified-fixed-workload") << "\n"
             << (routing_favored_workload
                 ? "  target mix: 20% sp-ok, 40% need-purify, 10% need-detour, 30% need-both; see selected mix above after fallback\n"
                 : "  mix (no purification): 50% robust, 30% mid-only, 20% low-only\n")
             << "  hop distribution:";
        for(const auto& [hop, count] : hop_dist)
            cerr << " " << hop << "hop=" << count;
        cerr << "\n================================================" << endl;
        assert((int)default_requests[r].size() == request_pool_size);
    }




    // vector<string> X_names = {"time_limit", "request_cnt", "num_nodes", "avg_memory", "tao"};
    vector<string> X_names = {"tao"};
    //vector<string> X_names = { "request_cnt", "time_limit", "tao",  "fidelity_threshold" , "avg_memory","hop_count","swap_prob" };
    //vector<string> X_names = {"Zmin","bucket_eps","time_eta"};
    vector<string> Y_names = {"fidelity_gain", "succ_request_cnt","actual_req_cnt"};
    vector<string> algo_names = {"ZFA_UB","ZFA2","ZFA_routing","MyAlgo1", "MyAlgo3","SP_skewed","SP_balanced"};
    // init result


    vector<PathMethod*> path_methods;
    path_methods.emplace_back(new Greedy());
    /* path_methods.emplace_back(new QCAST());
    path_methods.emplace_back(new REPS()); */
    for(PathMethod *path_method : path_methods) {

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
                  try {
                    cerr << "[CKPT] === ROUND " << r << " START | X=" << X_name << " val=" << change_value << " ===" << endl;
                    DBG_mem("round_start");
                    string filename = file_path + "input/round_" + to_string(r) + ".input";
                    ofstream ofs;
                    ofs.open(file_path + "log/" + path_method->get_name() + "_" + X_name + "_in_" + to_string(change_value) + "_Round_" + to_string(r) + ".log");

                    time_t now = time(0);
                    char* dt = ctime(&now);
                    cerr  << "時間 " << dt << endl << endl;
                    ofs << "時間 " << dt << endl << endl;




                    double A = 0.25, B = 0.75, tao = input_parameter["tao"], T = 10, n = 2;
                    DBG_HERE("before Graph ctor");
                    Graph graph(filename, time_limit, swap_prob, avg_memory, min_fidelity, max_fidelity, fidelity_threshold, A, B, n, T, tao,Zmin,bucket_eps,time_eta,input_parameter["delta_P"],input_parameter["entangle_lambda"],input_parameter["entangle_time"]);
                    DBG_HERE("after Graph ctor");
                    DBG_mem("after_graph");

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
                        DBG_HERE("requests filled from default_requests");
                    }
                    else{
                        DBG_HERE("before generate_requests_fid");
                        requests=generate_requests_fid(graph,request_cnt,0,hop_count);
                        DBG_HERE("after generate_requests_fid");
                    }
                    cerr << "[CKPT] requests.size()=" << requests.size() << endl;
                    DBG_HERE("before path_graph copy");
                    Graph path_graph = graph;
                    DBG_HERE("after path_graph copy");
                    DBG_mem("after_path_graph_copy");
                    path_graph.increase_resources(10);
                    DBG_HERE("after increase_resources");
                    PathMethod *new_path_method;
                    if(path_method->get_name() == "Greedy") new_path_method = new Greedy();
                    else if(path_method->get_name() == "QCAST") new_path_method = new QCAST();
                    else if(path_method->get_name() == "REPS") new_path_method = new REPS();
                    else {
                        cerr << "unknown path method" << endl;
                        assert(false);
                    }

                    DBG_HERE("before build_paths");
                    new_path_method->build_paths(path_graph, requests);
                    DBG_HERE("after build_paths");
                    DBG_mem("after_build_paths");
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
                    DBG_HERE("after path_st/paths build");

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
                    #ifndef WPFA_QUICK_BUILD
                    if(!compare_only) {
                        DBG_HERE("before new WernerAlgo3");
                        algorithms.emplace_back(new WernerAlgo3(graph,requests,paths));  // ZFA_UB (LP upper bound with purify)
                        DBG_HERE("after new WernerAlgo3");
                        DBG_mem("after_new_WernerAlgo3");
                    }
                    if(!compare_only) {
                        DBG_HERE("before new WernerAlgo2");
                        auto* zfa2 = new WernerAlgo2(graph,requests,paths);
                        DBG_HERE("after new WernerAlgo2");
                        DBG_mem("after_new_WernerAlgo2");
                        string exp_label = X_name + "=" + to_string(change_value) + " Round=" + to_string(r);
                        zfa2->set_experiment_label(exp_label);
                        algorithms.emplace_back(zfa2);
                    }
                    #endif
                    {
                        DBG_HERE("before new WernerAlgo_routing");
                        auto* zfa_routing = new WernerAlgo_routing(graph,requests,paths);
                        DBG_HERE("after new WernerAlgo_routing");
                        DBG_mem("after_new_WernerAlgo_routing");
                        string exp_label = X_name + "=" + to_string(change_value) + " Round=" + to_string(r);
                        zfa_routing->set_experiment_label(exp_label);
                        algorithms.emplace_back(zfa_routing);
                    }
                    if(X_name!="Zmin"&&X_name!="bucket_eps"&&X_name!="time_eta"){
                        #ifndef WPFA_QUICK_BUILD
                        if(!compare_only) {
                            DBG_HERE("before new MyAlgo1");
                            algorithms.emplace_back(new MyAlgo1(graph, requests, paths));
                            DBG_HERE("after new MyAlgo1");
                            DBG_HERE("before new MyAlgo3");
                            algorithms.emplace_back(new MyAlgo3(graph, requests, paths));
                            DBG_HERE("after new MyAlgo3");
                        }
                        #endif
                        // SP baselines：最短路徑 + 固定 swap 排程（skewed / balanced）
                        algorithms.emplace_back(new GreedySP(graph, requests, paths, GreedySP::SwapMode::SKEWED));
                        algorithms.emplace_back(new GreedySP(graph, requests, paths, GreedySP::SwapMode::BALANCED));
                        DBG_mem("after_all_algos_ctor");
                    }


                    // 統一設定實驗標籤（routing_trace.csv 的 exp_label 欄位）
                    {
                        string exp_label_all = X_name + "=" + to_string(change_value) + " Round=" + to_string(r);
                        for(auto* algo : algorithms)
                            algo->set_experiment_label(exp_label_all);
                    }

                    //#pragma omp parallel for schedule(dynamic)
                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        cerr << "[CKPT] >>> RUN algo[" << i << "] = " << algorithms[i]->get_name() << endl;
                        DBG_mem("before_run");
                        try {
                            algorithms[i]->run();
                        } catch(const std::bad_alloc& e) {
                            cerr << "[FATAL] std::bad_alloc inside algo[" << i << "] = "
                                 << algorithms[i]->get_name() << " : " << e.what() << endl;
                            DBG_mem("at_bad_alloc");
                            throw;
                        } catch(const std::exception& e) {
                            cerr << "[FATAL] std::exception inside algo[" << i << "] = "
                                 << algorithms[i]->get_name() << " : " << e.what() << endl;
                            throw;
                        }
                        cerr << "[CKPT] <<< DONE algo[" << i << "] = " << algorithms[i]->get_name() << endl;
                        DBG_mem("after_run");
                    }

                    // 該輪所有演算法跑完 → 按 request 分組寫出 routing trace（txt + csv）
                    AlgorithmBase::flush_routing_trace();



                    for(int i = 0; i < (int)algorithms.size(); i++) {
                        for(string Y_name : Y_names) {
                            result[r][algorithms[i]->get_name()][Y_name] = algorithms[i]->get_res(Y_name);
                        }
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
                    cerr << "[CKPT] === ROUND " << r << " END ===" << endl;
                    DBG_mem("round_end");
                  } catch(const std::bad_alloc& e) {
                      cerr << "[FATAL] std::bad_alloc in round r=" << r << " : " << e.what() << endl;
                      DBG_mem("round_bad_alloc");
                      throw;
                  } catch(const std::exception& e) {
                      cerr << "[FATAL] std::exception in round r=" << r << " : " << e.what() << endl;
                      throw;
                  }

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
            }
        }
    }
    return 0;
}
