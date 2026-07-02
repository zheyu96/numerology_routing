#ifndef __GREEDY_SP_H
#define __GREEDY_SP_H

#include "../AlgorithmBase/AlgorithmBase.h"
#include "../../Network/Graph/Graph.h"
#include "../../config.h"

#include <vector>
#include <map>
#include <utility>
#include <string>

using namespace std;

/**
 * GreedySP：最短路徑 + 固定 swap 排程的 greedy baseline
 *
 * 每條 request 依序（FCFS）處理：
 *   1. 在「未被 ban 的節點」子圖上 BFS 找最短路徑
 *   2. 在該路徑上用固定結構的 swap 排程蓋 connection：
 *        SKEWED   → 線性鏈（左深樹）：每個 timestep 往右延伸一段並 swap
 *        BALANCED → 平衡二元樹：兩半各自完成後在中點 swap
 *   3. purification 採全 link 均一輪數，從 0 輪開始試到 3 輪，
 *      取第一個過 fidelity threshold 的輪數（fidelity 模型與 ZFA 系列
 *      完全相同：Shape::get_fidelity + Purification pumping）
 *   4. 對完成時間 t_final 從最早可行掃到 time_limit-1，找第一個
 *      memory 足夠（含 purification 額外對數，Purify_in_vt 表與 ZFA 一致）
 *      的 slot；全部不行 → ban 掉 memory 瓶頸節點、換下一條路徑重試
 *   5. 重試 MAX_PATH_TRIES 次仍失敗 → 放棄該 request
 *
 * 資源預留與統計（reserve_shape / fidelity_gain / succ_request_cnt）
 * 與其他演算法走同一套 Graph 介面，保證比較公平。
 */
class GreedySP : public AlgorithmBase {
public:
    enum class SwapMode { SKEWED, BALANCED };

    GreedySP(const Graph& graph,
             const vector<pair<int,int>>& requests,
             const map<SDpair, vector<Path>>& paths,
             SwapMode mode);

    void run();

private:
    SwapMode swap_mode;

    static const int MAX_PURIFY_ROUNDS = 3;   // 與 ZFA 系列的 purify_time 一致
    static const int MAX_PATH_TRIES   = 5;    // 每條 request 最多換幾條路徑

    // 與 WernerAlgo 系列相同：purify r 輪時，link 兩端在第 i 個 offset
    // 需要的總 memory 數（-1 之後是額外量）
    double Purify_in_vt[4][5] = {
        {1,1},
        {1,2,2},
        {1,2,3,2},
        {1,2,3,3,2},
    };

    // 在排除 banned 節點的子圖上找 src->dst 的 hop 最短路徑（空 = 不連通）
    vector<int> bfs_path(int src, int dst, const vector<bool>& banned);

    // swap 樹的分割點與最早可行完成時間
    int split_point(int lo, int hi);
    int earliest_finish(int lo, int hi, const vector<int>& rounds);

    // 建構 links [lo,hi) 在時間 t 完成的 Shape_vector（格式同 backtrack_shape）
    Shape_vector build_segment(const vector<int>& path, int lo, int hi, int t,
                               const vector<int>& rounds);

    // 檢查 base + purification 額外 memory 是否足夠；不夠時回傳瓶頸節點
    bool memory_feasible(const Shape_vector& sv, const vector<int>& rounds,
                         int& bottleneck_node);

    // 預留 purification 額外的 memory（base 由 reserve_shape 扣）
    void reserve_purify_extra(const Shape_vector& sv, const vector<int>& rounds);
};

#endif // __GREEDY_SP_H
