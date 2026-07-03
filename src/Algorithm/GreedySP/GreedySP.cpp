#include "../../Network/Purification/Purification.h"
#include "GreedySP.h"
#include <iostream>
#include <queue>
#include <cmath>
#include <algorithm>

using namespace std;

GreedySP::GreedySP(const Graph& graph,
                   const vector<pair<int,int>>& requests,
                   const map<SDpair, vector<Path>>& paths,
                   SwapMode mode)
    : AlgorithmBase(graph, requests, paths), swap_mode(mode)
{
    algorithm_name = (mode == SwapMode::SKEWED) ? "SP_skewed" : "SP_balanced";
}

vector<int> GreedySP::bfs_path(int src, int dst, const vector<bool>& banned) {
    int V = graph.get_num_nodes();
    vector<int> parent(V, -1);
    vector<bool> vis(V, false);
    queue<int> que;
    vis[src] = true;
    que.push(src);
    while (!que.empty()) {
        int u = que.front(); que.pop();
        if (u == dst) break;
        for (int v : graph.adj_list[u]) {
            if (vis[v] || banned[v]) continue;
            vis[v] = true;
            parent[v] = u;
            que.push(v);
        }
    }
    if (!vis[dst]) return {};
    vector<int> path;
    for (int v = dst; v != -1; v = parent[v]) path.push_back(v);
    reverse(path.begin(), path.end());
    return path;
}

int GreedySP::split_point(int lo, int hi) {
    // SKEWED：左深樹，每次只切最右一條 link；BALANCED：對半切
    return (swap_mode == SwapMode::SKEWED) ? (hi - 1) : ((lo + hi) / 2);
}

int GreedySP::earliest_finish(int lo, int hi, const vector<int>& rounds) {
    // 葉子（單條 link）：entangle 區間 [t-r-1, t] 需 t-r-1 >= 1
    if (hi - lo == 1) return rounds[lo] + 2;
    int m = split_point(lo, hi);
    return max(earliest_finish(lo, m, rounds), earliest_finish(m, hi, rounds)) + 1;
}

Shape_vector GreedySP::build_segment(const vector<int>& path, int lo, int hi, int t,
                                     const vector<int>& rounds) {
    if (hi - lo == 1) {
        int r = rounds[lo];
        Shape_vector leaf;
        leaf.push_back({path[lo],     {{t - r - 1, t}}});
        leaf.push_back({path[lo + 1], {{t - r - 1, t}}});
        return leaf;
    }
    int m = split_point(lo, hi);
    Shape_vector left  = build_segment(path, lo, m, t - 1, rounds);
    Shape_vector right = build_segment(path, m, hi, t - 1, rounds);

    // 與 WernerAlgo backtrack_shape 的 MERGE 相同：junction 節點持有左右
    // 兩段的區間，最外側兩端在 swap 完成的時間 t 各延長一格
    Shape_vector result = left;
    result.back().second.push_back(right.front().second.front());
    for (int i = 1; i < (int)right.size(); i++) result.push_back(right[i]);
    result.front().second[0].second++;
    result.back().second[0].second++;
    return result;
}

bool GreedySP::memory_feasible(const Shape_vector& sv, const vector<int>& rounds,
                               int& bottleneck_node) {
    bottleneck_node = -1;
    map<pair<int,int>, int> total_need; // (node, t) -> amount

    // base：shape 每個 (node, timeslot) 各 1
    for (const auto& nd : sv) {
        for (const auto& rng : nd.second)
            for (int t = rng.first; t <= rng.second; t++)
                total_need[{nd.first, t}]++;
    }
    // purification 額外對數（offset 對應方式與 gen_leaf_label / rounding 一致）
    for (size_t li = 0; li + 1 < sv.size(); li++) {
        int r = (li < rounds.size()) ? rounds[li] : 0;
        if (r <= 0) continue;
        int link_start = sv[li].second.back().first;
        int u = sv[li].first, v = sv[li + 1].first;
        for (int ti = 0; ti <= r + 1; ti++) {
            int extra = (int)Purify_in_vt[r][r + 1 - ti] - 1;
            if (extra <= 0) continue;
            int t = link_start + ti;
            if (t >= time_limit) return false;
            total_need[{u, t}] += extra;
            total_need[{v, t}] += extra;
        }
    }
    for (const auto& [nt, amount] : total_need) {
        if (graph.get_node_memory_at(nt.first, nt.second) < amount) {
            bottleneck_node = nt.first;
            return false;
        }
    }
    return true;
}

void GreedySP::reserve_purify_extra(const Shape_vector& sv, const vector<int>& rounds) {
    for (size_t li = 0; li + 1 < sv.size(); li++) {
        int r = (li < rounds.size()) ? rounds[li] : 0;
        if (r <= 0) continue;
        int link_start = sv[li].second.back().first;
        int u = sv[li].first, v = sv[li + 1].first;
        for (int ti = 0; ti <= r + 1; ti++) {
            int extra = (int)Purify_in_vt[r][r + 1 - ti] - 1;
            if (extra <= 0) continue;
            int t = link_start + ti;
            graph.reserve_node_memory_at(u, t, extra);
            graph.reserve_node_memory_at(v, t, extra);
        }
    }
}

void GreedySP::run() {
    int V = graph.get_num_nodes();
    double fid_th = graph.get_fidelity_threshold();
    vector<int> finished;
    int fail_fid = 0, fail_mem = 0;

    for (int i = 0; i < (int)requests.size(); i++) {
        int src = requests[i].first, dst = requests[i].second;
        vector<bool> banned(V, false);
        bool served = false;
        bool fid_hopeless = false;

        for (int attempt = 0; attempt < MAX_PATH_TRIES && !served; attempt++) {
            vector<int> path = bfs_path(src, dst, banned);
            if (path.empty()) break;
            int h = (int)path.size() - 1;

            // 1) 找最小的均一 purify 輪數 r 使 fidelity 過 threshold。
            //    fidelity 只跟結構與時間差有關（time-shift invariant），
            //    用最早完成時間探測一次即可。
            int chosen_r = -1;
            for (int r = 0; r <= MAX_PURIFY_ROUNDS && chosen_r < 0; r++) {
                vector<int> rounds(h, r);
                int min_t = earliest_finish(0, h, rounds);
                if (min_t > time_limit - 1) continue; // 這個 r 連時間都不夠
                Shape_vector sv = build_segment(path, 0, h, min_t, rounds);
                double fid = 0;
                try {
                    Shape probe(sv, rounds);
                    fid = probe.get_fidelity(A, B, n, T, tao, graph.get_F_init(), true);
                } catch (const runtime_error&) {
                    fid = 0;
                }
                if (fid + EPS >= fid_th) chosen_r = r;
            }
            if (chosen_r < 0) {
                // 最短路都過不了 threshold，更長的替代路徑只會更差 → 放棄
                fid_hopeless = true;
                break;
            }

            // 2) 掃完成時間找第一個 memory 塞得下的 slot
            vector<int> rounds(h, chosen_r);
            int min_t = earliest_finish(0, h, rounds);
            int bottleneck = -1;
            for (int tf = min_t; tf <= time_limit - 1 && !served; tf++) {
                Shape_vector sv = build_segment(path, 0, h, tf, rounds);
                int bn = -1;
                if (!memory_feasible(sv, rounds, bn)) { bottleneck = bn; continue; }

                Shape shape(sv, rounds);
                graph.reserve_shape(shape, true);
                reserve_purify_extra(sv, rounds);
                served = true;

                // routing trace CSV（per-request 對照用）
                {
                    double fid_acc = shape.get_fidelity(A, B, n, T, tao, graph.get_F_init(), true);
                    double prob_acc = graph.path_Pr_purify(shape);
                    log_routing_trace(i, src, dst, "accepted", sv, rounds, fid_acc, prob_acc);
                }

                cerr << "\033[1;36m" << "[" << algorithm_name << " path] req#" << i
                     << " (" << src << "->" << dst << "): ";
                for (int j = 0; j <= h; j++) {
                    cerr << path[j];
                    if (j < h) cerr << " --(pur=" << chosen_r << ")--> ";
                }
                cerr << " | hop=" << h << " finish_t=" << tf << "\033[0m" << endl;
            }

            // 3) memory 不夠 → ban 掉瓶頸節點換路；瓶頸是端點就沒救了
            if (!served) {
                if (bottleneck == -1 || bottleneck == src || bottleneck == dst) break;
                banned[bottleneck] = true;
            }
        }

        if (served) finished.push_back(i);
        else if (fid_hopeless) fail_fid++;
        else fail_mem++; // 含找不到路徑的情況

        if (!served) {
            log_routing_trace(i, src, dst, fid_hopeless ? "fail_fid" : "fail_mem",
                              Shape_vector{}, vector<int>{}, 0, 0);
        }
    }

    cerr << "\033[1;35m" << "[" << algorithm_name << "] served=" << finished.size()
         << "/" << requests.size()
         << " | fail_fidelity=" << fail_fid
         << " | fail_mem_or_path=" << fail_mem
         << "\033[0m" << endl;

    sort(finished.rbegin(), finished.rend());
    for (int fin : finished) requests.erase(requests.begin() + fin);

    update_res();
    cerr << "[" << algorithm_name << "] end" << endl;
}
