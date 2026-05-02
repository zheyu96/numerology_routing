#include "Graph.h"
#include "../Purification/Purification.h"
#include <cmath>
using namespace std;

int rnd(int lower_bound, int upper_bound) {
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);
    return unif(generator) % (upper_bound - lower_bound + 1) + lower_bound;
}

Graph::Graph(string filename, int _time_limit, double _swap_prob, int avg_memory, double min_fidelity, double max_fidelity, double _fidelity_threshold, double _A, double _B, double _n, double _T, double _tao,double _Zmin,double _bucket_eps,double _time_eta,double _delta_P):
    time_limit(_time_limit), fidelity_threshold(_fidelity_threshold), A(_A), B(_B), n(_n), T(_T), tao(_tao),Zmin(_Zmin),bucket_eps(_bucket_eps),time_eta(_time_eta), fidelity_gain(0), usage(0), succ_request_cnt(0), actual_succ_request_cnt(0), delta_P(_delta_P) {
    // geneator an adj list

    ifstream graph_file(filename);
    if (!graph_file) {
        throw runtime_error("Cannot open input file: " + filename);
    }

    file_name = filename;
    graph_file >> num_nodes;
    
    adj_set.clear();
    adj_set.resize(num_nodes);
    adj_list.clear();
    adj_list.resize(num_nodes);

    swapping_succ_prob = _swap_prob;
    for(int id = 0; id < num_nodes; id++) {
        int memory_rand;
        double swap_prob = swapping_succ_prob;
        graph_file >> memory_rand;
        memory_rand += avg_memory;
        nodes.push_back(Node(id, memory_rand, time_limit, swap_prob));
    }
    int num_edges;
    graph_file >> num_edges;
    avg_entangle_prob = 0;
    for(int i = 0; i < num_edges; i++) {
        int v, u;
        double f_init,entangle_prob;
        graph_file >> v >> u >> entangle_prob;
        f_init=entangle_prob;
        assert(v != u);
        f_init = (f_init) * (max_fidelity - min_fidelity) + min_fidelity;
        adj_list[v].push_back(u);
        adj_list[u].push_back(v);
        F_init[{v, u}] = f_init;
        F_init[{u, v}] = f_init;
        entangle_succ_prob[{v, u}] = entangle_prob;
        entangle_succ_prob[{u, v}] = entangle_prob;
        avg_entangle_prob += entangle_prob;
    }
    avg_entangle_prob /= num_edges;

    for(int i = 0; i < num_nodes; i++) {
        for(auto v : adj_list[i]) {
            adj_set[i].insert(v);
        }
    }

    boundary.clear();
    for(double d = 0.9; d <= 1 + EPS; d += 0.01) {
        boundary.push_back(d);
    }

    cnt.clear();
    cnt.resize(boundary.size(), 0);
}

Graph::~Graph() { if(DEBUG) cerr << "Delete Graph" << endl; }

int Graph::get_node_memory_at(int node_id, int t) { return nodes[node_id].get_memory_at(t); }
int Graph::get_node_memory(int node_id) { return nodes[node_id].get_memory(); }
double Graph::get_node_swap_prob(int node_id) { return nodes[node_id].get_swap_prob(); }
int Graph::get_num_nodes() { return num_nodes; }
int Graph::get_time_limit() { return time_limit; }

double Graph::get_A() { return A; }
double Graph::get_B() { return B; }
double Graph::get_n() { return n; }
double Graph::get_T() { return T; }
double Graph::get_tao() { return tao; }
double Graph::get_entangle_succ_prob(int u, int v) { return entangle_succ_prob[{u, v}]; };
double Graph::get_fidelity_gain() { return fidelity_gain; }
double Graph::get_pure_fidelity() { return pure_fidelity; }
double Graph::get_fidelity_threshold() { return fidelity_threshold; }
double Graph::get_Zmin(){ return Zmin; }
double Graph::get_bucket_eps(){ return bucket_eps; }
double Graph::get_time_eta(){ return time_eta; }
double Graph::get_succ_request_cnt() { return succ_request_cnt;}
int Graph::get_actual_succ_request_cnt() { return actual_succ_request_cnt; }
int Graph::get_usage() { return usage; }

const vector<double>& Graph::get_boundary() const { return boundary; }
const vector<double>& Graph::get_cnt() const { return cnt; }

double Graph::get_F_init(int u, int v) {
    assert(adj_set[u].count(v));
    return F_init[{u, v}];
}

const map<pair<int, int>, double>& Graph::get_F_init() const { return F_init; }

double Graph::get_link_werner(int u,int v){ //werner state
    auto it=F_init.find({u,v});
    assert(it!=F_init.end());
    double W=Purification::fidelity_to_werner(it->second);
    assert(W>=0.0&&W<=1.0);
    return W;
}
double Graph::get_edge_W(int u,int v){ //werner log
    double w=get_link_werner(u,v);
    double eps=1e-12;
    return -log(max(w,eps));
}
void DFS(int x, vector<bool> &vis, vector<int> &par, vector<vector<int>> &adj) {
    vis[x] = true;
    for(auto v : adj[x]) {
        if(!vis[v]) {
            par[v] = x;
            DFS(v, vis, par, adj);
        }
    }
}

Path Graph::get_path(int from, int to) {
    vector<bool> vis(num_nodes + 1, false);
    vector<int> par(num_nodes + 1, -1);
    par[from] = -1;
    DFS(from, vis, par, adj_list);
    vector<int> path;
    int cur = to;
    while(cur != -1) {
        path.push_back(cur);
        cur = par[cur];
    }
    reverse(path.begin(), path.end());

    return path;
}

int Graph::distance(int src, int dst) {
    vector<bool> vis(num_nodes + 1, false);
    vector<int> dis(num_nodes + 1, 2000);
    // BFS
    queue<int> que;
    dis[src] = 0;
    vis[src] = true;
    que.push(src);
    while(!que.empty()) {
        int frt = que.front();
        que.pop();
        for(int v : adj_list[frt]) {
            if(!vis[v]) {
                que.push(v);
                dis[v] = min(dis[v], dis[frt] + 1);
                if(v == dst) return dis[v];
            }
        }
    }
    assert(false);
    return 2000;
}
double Graph::get_ini_fid(int src,int dst){
    vector<bool> vis(num_nodes+1,false);
    vector<double> f(num_nodes+1,0.0);
    f[src]=1.0;
    queue<int> que;
    que.push(src);
    while(que.size()){
        int cur=que.front();
        que.pop();
        if(vis[cur]) continue;
        vis[cur]=1;
        for(int v:adj_list[cur]){
            if(!vis[v]){
                double  fid=f[cur]*get_F_init(cur,v);
                if(fid>f[v]){
                    f[v]=fid;
                    que.push(v);
                }
            }
        }
    }
    return f[dst];
}
double Graph::get_delta_P(){
    return delta_P;    
}
// [修改] 增加 enable_purification 參數
bool Graph::check_resource(Shape shape, bool threshold /*= true*/, bool enable_purification /*= false*/) {
    Shape_vector nm = shape.get_node_mem_range();
    // [修改] 將參數傳遞給 get_fidelity
    if(threshold && shape.get_fidelity(A, B, n, T, tao, F_init, enable_purification) < fidelity_threshold) return false;
    for(int i = 0; i < (int)nm.size(); i++) {
        int node = nm[i].first;
        map<int, int> need_amount; // time to amount
        for(pair<int, int> rng : nm[i].second) {
            int left = rng.first, right = rng.second;
            if(left < 0) {
                cerr << "the reserve time is negtive" << endl;
                exit(1);
            }
            if(right >= time_limit) {
                cerr << "the reserve time is exceed the timelimit" << endl;
                cerr << "timelimt = " << time_limit << " reserve time = " << right << endl;
                exit(1);
            }
            for(int t = left; t <= right; t++) {
                need_amount[t]++;
            }
        }
        for(auto P : need_amount) {
            int t = P.first, amount = P.second;
            if(nodes[node].get_memory_at(t) < amount) {
                return false;
            }
        }
    }
    return true;
}

// [修改] 增加 enable_purification 參數 (雖然目前 ASAP 演算法可能用不到，但保持一致性)
bool Graph::check_resource_ASAP(Shape shape, bool threshold, bool enable_purification) {
    Shape_vector nm = shape.get_node_mem_range();
    if(threshold && shape.get_fidelity(A, B, n, T, tao, F_init, enable_purification) < fidelity_threshold) return false;
    int mx_amount = 0, earliest = time_limit, lastest = 0; 
    for(int i = 0; i < (int)nm.size(); i++) {
        map<int, int> need_amount; // time to amount
        for(pair<int, int> rng : nm[i].second) {
            int left = rng.first, right = rng.second;
            if(left < 0) {
                cerr << "the reserve time is negtive" << endl;
                exit(1);
            }
            if(right >= time_limit) {
                cerr << "the reserve time is exceed the timelimit" << endl;
                cerr << "timelimt = " << time_limit << " reserve time = " << right << endl;
                exit(1);
            }
            for(int t = left; t <= right; t++) {
                need_amount[t]++;
            }
        }
        for(auto P : need_amount) {
            int t = P.first, amount = P.second;
            earliest = min(earliest, t);
            lastest = max(lastest, t);
            mx_amount = max(mx_amount, amount);
        }
    }
    for(int i = 0; i < (int)nm.size(); i++) {
        int node = nm[i].first;
        for(int t = earliest; t <= lastest; t++) {
            if(nodes[node].get_memory_at(t) < mx_amount) {
                return false;
            }
        }
    }
    return true;
}

// [修改] 增加 enable_purification 參數
void Graph::reserve_shape_ASAP(Shape shape, bool enable_purification) {
    shape.check_valid();
    // cerr << "checked" << endl;
    Shape_vector nm = shape.get_node_mem_range();
    int mx_amount = 0, earliest = time_limit, lastest = 0; 
    for(int i = 0; i < (int)nm.size(); i++) {
        map<int, int> need_amount; // time to amount
        for(pair<int, int> rng : nm[i].second) {
            int left = rng.first, right = rng.second;

            if(right >= time_limit) {
                cerr << "the reserve time is exceed the timelimit" << endl;
                cerr << "timelimt = " << time_limit << " reserve time = " << right << endl;
                assert(false);
                exit(1);
            }
            for(int t = left; t <= right; t++) {
                need_amount[t]++;
            }
        }
        for(auto P : need_amount) {
            int t = P.first, amount = P.second;
            earliest = min(earliest, t);
            lastest = max(lastest, t);
            mx_amount = max(mx_amount, amount);
        }
    }
    for(int i = 0; i < (int)nm.size(); i++) {
        int node = nm[i].first;
        for(int t = earliest; t <= lastest; t++) {
            if(nodes[node].get_memory_at(t) < mx_amount) {
                cerr << "node " << node << "\'s memory is not enough at time " << t << endl;
                exit(1);
            }
            usage += mx_amount;
            nodes[node].reserve_memory(t, mx_amount);
        }
    }
    for(int i = 1; i < (int)nm.size(); i++) {
        int node1 = nm[i - 1].first;
        int node2 = nm[i].first;
        if(adj_set[node1].count(node2) == 0) {
            cerr << "shape error, the next node is not connected" << endl;
            exit(1);
        }
    }

    double shape_fidelity = shape.get_fidelity(A, B, n, T, tao, F_init, enable_purification);
    if(shape_fidelity > fidelity_threshold) {
        int src=nm.front().first,dst=nm.back().first;
        double pr = enable_purification ? path_Pr_purify(shape) : path_Pr(shape);
        cerr<<"\033[1;34m"<<"initial fid = "<<get_ini_fid(src,dst)<<" real fidelity = "<<shape_fidelity<<" hop count = "<<(nm.size()-1)<<" prob = "<<pr<<"\033[0m"<<endl;
        fidelity_gain += (Purification::fidelity_to_werner(shape_fidelity) * pr);
        pure_fidelity += shape_fidelity;
        succ_request_cnt += pr;
        actual_succ_request_cnt++;
    }

    for(int i = 0; i < (int)boundary.size(); i++) {
        if(shape_fidelity < boundary[i]) {
            cnt[i] = cnt[i] + 1;
            break;
        }
    }
}

// [修改] 增加 enable_purification 參數
void Graph::reserve_shape(Shape shape, bool enable_purification /*= false*/) {
    shape.check_valid();
    // cerr << "checked" << endl;
    Shape_vector nm = shape.get_node_mem_range();
    for(int i = 0; i < (int)nm.size(); i++) {
        int node = nm[i].first;
        map<int, int> need_amount; // time to amount
        for(pair<int, int> rng : nm[i].second) {
            int left = rng.first, right = rng.second;

            if(right >= time_limit) {
                cerr << "the reserve time is exceed the timelimit" << endl;
                cerr << "timelimt = " << time_limit << " reserve time = " << right << endl;
                assert(false);
                exit(1);
            }
            for(int t = left; t <= right; t++) {
                need_amount[t]++;
            }
        }
        for(auto P : need_amount) {
            int t = P.first, amount = P.second;
            if(nodes[node].get_memory_at(t) < amount) {
                cerr << "node " << node << "\'s memory is not enough at time " << t << endl;
                exit(1);
            }
            usage += amount;
            nodes[node].reserve_memory(t, amount);
        }
    }
    for(int i = 1; i < (int)nm.size(); i++) {
        int node1 = nm[i - 1].first;
        int node2 = nm[i].first;
        if(adj_set[node1].count(node2) == 0) {
            cerr << "shape error, the next node is not connected" << endl;
            exit(1);
        } 
    }

    // [修改] 計算 fidelity 時使用 enable_purification
    double shape_fidelity = shape.get_fidelity(A, B, n, T, tao, F_init, enable_purification);
    
    if(shape_fidelity + EPS < fidelity_threshold) {
        cerr << "the fidelity of shape is not greater than threshold" << endl;
        // Debug 訊息
        cerr << "Shape Fid: " << shape_fidelity << ", Thres: " << fidelity_threshold << endl;
        assert(false);
        exit(1);
    }
    int src=nm.front().first,dst=nm.back().first;
    double pr = enable_purification ? path_Pr_purify(shape) : path_Pr(shape);
    cerr<<"\033[1;34m"<<"initial fid = "<<get_ini_fid(src,dst)<<" real fidelity = "<<shape_fidelity<<" hop count = "<<(nm.size()-1)<<" prob = "<<pr<<"\033[0m"<<endl;
    fidelity_gain += (Purification::fidelity_to_werner(shape_fidelity) * pr);
    pure_fidelity += shape_fidelity;
    succ_request_cnt += pr;
    actual_succ_request_cnt++;

    for(int i = 0; i < (int)boundary.size(); i++) {
        if(shape_fidelity < boundary[i]) {
            cnt[i] = cnt[i] + 1;
            break;
        }
    }
}

// [修改] 增加 enable_purification 參數
void Graph::reserve_shape2(Shape shape, bool enable_purification) {
    shape.check_valid();
    // cerr << "checked" << endl;
    Shape_vector nm = shape.get_node_mem_range();
    for(int i = 0; i < (int)nm.size(); i++) {
        int node = nm[i].first;
        map<int, int> need_amount; // time to amount
        for(pair<int, int> rng : nm[i].second) {
            int left = rng.first, right = rng.second;

            if(right >= time_limit) {
                cerr << "the reserve time is exceed the timelimit" << endl;
                cerr << "timelimt = " << time_limit << " reserve time = " << right << endl;
                assert(false);
                exit(1);
            }
            for(int t = left; t <= right; t++) {
                need_amount[t]++;
            }
        }
        for(auto P : need_amount) {
            int t = P.first, amount = P.second;
            if(nodes[node].get_memory_at(t) < amount) {
                cerr << "node " << node << "\'s memory is not enough at time " << t << endl;
                exit(1);
            }
            usage += amount;
            nodes[node].reserve_memory(t, amount);
        }
    }
    for(int i = 1; i < (int)nm.size(); i++) {
        int node1 = nm[i - 1].first;
        int node2 = nm[i].first;
        if(adj_set[node1].count(node2) == 0) {
            cerr << "shape error, the next node is not connected" << endl;
            exit(1);
        } 
    }

    double shape_fidelity = shape.get_fidelity(A, B, n, T, tao, F_init, enable_purification);
    if(shape_fidelity > fidelity_threshold) {
        int src=nm.front().first,dst=nm.back().first;
        double pr = enable_purification ? path_Pr_purify(shape) : path_Pr(shape);
        cerr<<"\033[1;34m"<<"initial fid = "<<get_ini_fid(src,dst)<<" real fidelity = "<<shape_fidelity<<" hop count = "<<(nm.size()-1)<<" prob = "<<pr<<"\033[0m"<<endl;
        fidelity_gain += (Purification::fidelity_to_werner(shape_fidelity) * pr);
        pure_fidelity += shape_fidelity;
        succ_request_cnt += pr;
        actual_succ_request_cnt++;
    }

    for(int i = 0; i < (int)boundary.size(); i++) {
        if(shape_fidelity < boundary[i]) {
            cnt[i] = cnt[i] + 1;
            break;
        }
    }
}
double Graph::path_Pr(Path path) {
    double Pr = 1;
    for(int i = 1; i < (int)path.size() - 1; i++) {
        Pr *= nodes[path[i]].get_swap_prob();
    }
    for(int i = 1; i < (int)path.size(); i++) {
        int node1 = path[i - 1], node2 = path[i];
        Pr *= get_entangle_succ_prob(node1, node2);
    }
    return Pr;
}
double Graph::path_Pr(const Shape& shape) {
    Path path;
    for(const auto& P : shape.get_node_mem_range()) {
        path.push_back(P.first);
    }
    return path_Pr(path);
}

double Graph::path_Pr_purify(const Shape& shape) {
    const Shape_vector& nm = shape.get_node_mem_range();
    const vector<int>& pr = shape.get_link_purify_rounds();
    double Pr = 1.0;

    // 每條 link 的 entanglement + purification 成功機率
    for(int i = 0; i < (int)nm.size() - 1; i++) {
        int u = nm[i].first, v = nm[i + 1].first;
        int rounds = (i < (int)pr.size()) ? pr[i] : 0;

        double fresh_werner = get_link_werner(u, v);
        double p_link = Purification::pumping_success_prob(
            get_entangle_succ_prob(u, v), fresh_werner, rounds);
        Pr *= p_link;
    }

    // 中間節點的 swap 成功機率
    for(int i = 1; i < (int)nm.size() - 1; i++) {
        Pr *= nodes[nm[i].first].get_swap_prob();
    }

    return Pr;
}

void Graph::reserve_path(Path path, int amount) {
    int src = path[0], dst = path.back();
    nodes[src].reserve_memory(-amount);
    nodes[dst].reserve_memory(-amount);
    for(int node : path) {
        if(nodes[node].get_memory() < 2 * amount) {
            cerr << "error: Node Memeory is not enough in reserve_path" << endl;
            cerr << "Node Memory = " << nodes[node].get_memory() << " amount = " << amount << endl;
            exit(1);
        }
        nodes[node].reserve_memory(2 * amount);
    }
}
void Graph::reserve_path(Path path) {
    int src = path[0], dst = path.back();
    int amount = min(nodes[src].get_memory(), nodes[dst].get_memory());
    for(int node : path) {
        if(node == src || node == dst) continue;
        amount = min(amount, nodes[node].get_memory() / 2);
    }

    if(amount == 0) {
        cerr << "error: in reserve_path memory is not enough" << endl;
        exit(1);
    }

    amount = 1;
    nodes[src].reserve_memory(amount);
    nodes[dst].reserve_memory(amount);
    for(int node : path) {
        if(node == src || node == dst) continue;
        nodes[node].reserve_memory(amount * 2);
    }
}
bool Graph::check_path_resource(Path path, int amount) {
    int src = path[0], dst = path.back();
    int remain = min(nodes[src].get_memory(), nodes[dst].get_memory());
    for(int node : path) {
        if(node == src || node == dst) continue;
        remain = min(remain, nodes[node].get_memory() / 2);
    }
    return remain >= amount;
}


void Graph::reserve_node_memory_at(int node_id, int t, int amount) {
    nodes[node_id].reserve_memory(t, amount);
    usage += amount;
}

void Graph::increase_resources(int multi) {
    for(int i = 0; i < num_nodes; i++) {
        int node_memory = nodes[i].get_memory();
        nodes[i].reserve_memory(-node_memory * (multi - 1));
    }
}
