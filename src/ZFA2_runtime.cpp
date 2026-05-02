#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <random>
#include <algorithm>
#include <chrono>
#include <omp.h>
#include <cassert>

#include "./config.h"
#include "Network/Graph/Graph.h"
#include "Algorithm/AlgorithmBase/AlgorithmBase.h"
#include "Algorithm/WernerAlgo2/WernerAlgo2.h"
#include "Network/PathMethod/Greedy/Greedy.h"

using namespace std;

SDpair generate_new_request(int num_of_node) {
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, num_of_node - 1);
    int node1 = unif(generator), node2 = unif(generator);
    while(node1 == node2) node2 = unif(generator);
    return make_pair(node1, node2);
}

vector<SDpair> generate_requests_fid(Graph graph, int requests_cnt, double fid_th, double hop_th) {
    int n = graph.get_num_nodes();
    vector<pair<SDpair, double>> cand[22];
    random_device rd;
    default_random_engine generator = default_random_engine(rd());
    uniform_int_distribution<int> unif(0, 1e9);

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i == j) continue;
            double fid = graph.get_ini_fid(i, j);
            assert(fid >= 0.0 && fid <= 1.0);
            if(fid > fid_th && graph.distance(i, j) >= hop_th) {
                int index = fid / 0.05;
                if(index < 0) continue;
                if(index > 20) index = 20;
                cand[index].emplace_back(make_pair(make_pair(i, j), (double)graph.distance(i, j)));
            }
        }
    }

    for(int i = 0; i < 22; i++) {
        shuffle(cand[i].begin(), cand[i].end(), generator);
    }

    vector<SDpair> requests;
    int pos[22] = {0};
    int idx = 0;
    while(requests.size() < (size_t)requests_cnt) {
        int cnt = unif(generator) % 5 + 4;
        cnt = min(cnt, (int)(requests_cnt - requests.size()));
        int try_idx = 0;
        while(cand[21 - idx].empty() && try_idx < 22) {
            idx = (idx + 1) % 22;
            try_idx++;
        }
        if(try_idx >= 22) break;

        for(int i = 0; i < cnt; i++) {
            requests.push_back(cand[21 - idx][pos[21 - idx]].first);
            pos[21 - idx] = (pos[21 - idx] + 1) % cand[21 - idx].size();
        }
        idx = (idx + 1) % 22;
    }
    return requests;
}

bool is_better_result(const map<string, double>& lhs, const map<string, double>& rhs) {
    const double EPS = 1e-9;

    double lhs_fidelity = lhs.at("fidelity_gain");
    double rhs_fidelity = rhs.at("fidelity_gain");
    if(lhs_fidelity > rhs_fidelity + EPS) return true;
    if(rhs_fidelity > lhs_fidelity + EPS) return false;

    double lhs_succ = lhs.at("succ_request_cnt");
    double rhs_succ = rhs.at("succ_request_cnt");
    if(lhs_succ > rhs_succ + EPS) return true;
    if(rhs_succ > lhs_succ + EPS) return false;

    return lhs.at("runtime") + EPS < rhs.at("runtime");
}

int main() {
    string file_path = "../data/";

    system("mkdir -p ../data/ans");
    system("mkdir -p ../data/log");
    system("mkdir -p ../data/input");
    {
        ofstream ofs(file_path + "ans/Greedy_ZFA2_actual_succ_request_cnt.txt", ios::out);
        ofs << "X_name value actual_succ_request_cnt" << endl;
    }

    map<string, double> default_setting;
    default_setting["num_nodes"] = 100;
    default_setting["request_cnt"] = 300;
    default_setting["time_limit"] = 13;
    default_setting["avg_memory"] = 6;
    default_setting["tao"] = 0.002;
    default_setting["fidelity_threshold"] = 0.7;
    default_setting["epsilon"] = 0.35;
    default_setting["Zmin"] = 0.02702867239;
    default_setting["bucket_eps"] = 0.01;
    default_setting["time_eta"] = 0.001;
    default_setting["delta_P"] = 0.01;

    map<string, vector<double>> change_parameter;
    change_parameter["epsilon"] = {0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4, 0.45,0.5, 0.55,0.6, 0.65,0.7, 0.75,0.8, 0.85,0.9, 0.95};
    change_parameter["bucket_eps"] = {0.001, 0.01,0.02,0.03,0.04, 0.05, 0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45, 0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95};

    vector<string> X_names = {"epsilon", "bucket_eps"};
    int round = 5;

    vector<vector<SDpair>> default_requests(round);
    cerr << "Pre-generating requests..." << endl;
    for(int r = 0; r < round; r++) {
        string filename = file_path + "input/round_" + to_string(r) + ".input";

        ifstream fin(filename);
        if(!fin.is_open()) {
            cerr << "Warning: Input file " << filename << " not found! Running generator..." << endl;
            string cmd = "python3 graph_generator.py " + filename + " " + to_string((int)default_setting["num_nodes"]);
            system(cmd.c_str());
        }

        Graph temp_g(filename, 20, 0.9, 10, 0.5, 0.98, 0.7, 0.25, 0.75, 2, 10,
                     0.002, 0.027, 0.01, 0.001, 0.01);
        default_requests[r] = generate_requests_fid(temp_g, (int)default_setting["request_cnt"], 0.7, 3);
        assert(!default_requests[r].empty());
    }

    for(const string& X_name : X_names) {
        if(change_parameter.find(X_name) == change_parameter.end()) continue;

        for(double change_value : change_parameter[X_name]) {
            cerr << "Testing " << X_name << " (WernerAlgo2 runtime) = " << change_value << endl;
            vector<map<string, map<string, double>>> result(round);

            #pragma omp parallel for
            for(int r = 0; r < round; r++) {
                string filename = file_path + "input/round_" + to_string(r) + ".input";

                map<string, double> input_param = default_setting;
                input_param[X_name] = change_value;

                Graph graph(filename,
                            (int)input_param["time_limit"], 0.9,
                            (int)input_param["avg_memory"], 0.5, 0.98,
                            input_param["fidelity_threshold"], 0.25, 0.75, 2, 10,
                            input_param["tao"], input_param["Zmin"],
                            input_param["bucket_eps"], input_param["time_eta"],
                            input_param["delta_P"]);

                vector<pair<int, int>> requests = default_requests[r];

                Greedy path_method;
                path_method.build_paths(graph, requests);
                const auto& paths = path_method.get_paths();

                AlgorithmBase* algo = new WernerAlgo2(graph, requests, paths,
                                                      input_param["epsilon"],
                                                      input_param["bucket_eps"]);

                auto start_time = chrono::high_resolution_clock::now();
                algo->run();
                auto end_time = chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed = end_time - start_time;

                string algo_name = algo->get_name();
                #pragma omp critical
                {
                    result[r][algo_name]["fidelity_gain"] = algo->get_res("fidelity_gain");
                    result[r][algo_name]["succ_request_cnt"] = algo->get_res("succ_request_cnt");
                    result[r][algo_name]["actual_succ_request_cnt"] = algo->get_res("actual_succ_request_cnt");
                    result[r][algo_name]["runtime"] = elapsed.count();
                }
                delete algo;
            }

            string target_algo = "ZFA2";
            int best_round = 0;
            for(int r = 1; r < round; r++) {
                if(is_better_result(result[r][target_algo], result[best_round][target_algo])) {
                    best_round = r;
                }
            }

            string ans_filename = "ans/Greedy_" + X_name + "_ZFA2.ans";
            ofstream ofs(file_path + ans_filename, ios::app);
            ofs << change_value << " "
                << result[best_round][target_algo]["fidelity_gain"] << " "
                << result[best_round][target_algo]["succ_request_cnt"] << " "
                << result[best_round][target_algo]["runtime"] << endl;
            ofs.close();

            ofstream actual_ofs(file_path + "ans/Greedy_ZFA2_actual_succ_request_cnt.txt", ios::app);
            actual_ofs << X_name << " "
                       << change_value << " "
                       << result[best_round][target_algo]["actual_succ_request_cnt"] << endl;
            actual_ofs.close();
        }
    }
    return 0;
}
