# --- 路徑與變數設定 ---
# 基礎依賴物件
dependency = Network/Node/Node.o Network/Shape/Shape.o Network/Graph/Graph.o Algorithm/AlgorithmBase/AlgorithmBase.o config.o

# 完整演算法物件 (用於 a.out / main.cpp)
algorithm = Algorithm/MyAlgo1/MyAlgo1.o Algorithm/MyAlgo2/MyAlgo2.o Algorithm/MyAlgo3/MyAlgo3.o \
            Algorithm/MyAlgo4/MyAlgo4.o Algorithm/MyAlgo5/MyAlgo5.o Algorithm/MyAlgo6/MyAlgo6.o \
            Algorithm/WernerAlgo/WernerAlgo.o Algorithm/WernerAlgo2/WernerAlgo2.o \
            Algorithm/WernerAlgo_routing/WernerAlgo_routing.o \
            Algorithm/WernerAlgo3/WernerAlgo3.o \
            Algorithm/WernerAlgo_UB/WernerAlgo_UB.o \
            Algorithm/GreedySP/GreedySP.o

# 時間量測專用演算法物件
algorithm_zfa = Algorithm/WernerAlgo_time/WernerAlgo_time.o 
algorithm_zfa2 = Algorithm/WernerAlgo2/WernerAlgo2.o
algorithm2 = Algorithm/MyAlgo3/MyAlgo3.o Algorithm/MyAlgo1_time/MyAlgo1_time.o

# 路由與路徑方法
routing = Network/PathMethod/PathMethodBase/PathMethod.o Network/PathMethod/Greedy/Greedy.o \
          Network/PathMethod/QCAST/QCAST.o Network/PathMethod/REPS/REPS.o

# 編譯旗標
gcc_flag = -O3 -fopenmp -march=native --std=c++17

# Gurobi 設定 (如果沒有安裝 Gurobi，留空即可，REPS 會編譯失敗但不影響主程式)
# gurobi_path = /home/oba/danger/gurobi951/linux64
# gurobi_flag = -I$(gurobi_path)/include/ -L$(gurobi_path)/lib/ -lgurobi_c++ -lgurobi95
gurobi_flag =

# --- 主要執行目標 ---

all: a.out

# 1. 執行 WernerAlgo_time 測試 (ZFA_runtime.cpp)
zfa_run: ZFA_runtime.cpp $(dependency) $(algorithm_zfa) $(routing)
	g++ -fopenmp -O3 $(gcc_flag) ZFA_runtime.cpp $(dependency) $(algorithm_zfa) $(routing) $(gurobi_flag) -o zfa_run
	./zfa_run 

# 2. 執行 WernerAlgo2 runtime 測試 (ZFA2_runtime.cpp)
zfa2_run: ZFA2_runtime.cpp $(dependency) $(algorithm_zfa2) $(routing)
	g++ -fopenmp -O3 $(gcc_flag) ZFA2_runtime.cpp $(dependency) $(algorithm_zfa2) $(routing) $(gurobi_flag) -o zfa2_run
	./zfa2_run

# 3. 舊有的時間量測目標
runtime: runtime.cpp $(dependency) $(algorithm2) $(routing)
	g++ -fopenmp $(gcc_flag) runtime.cpp $(dependency) $(algorithm2) $(routing) $(gurobi_flag) -o runtime

# 4. 預設主程式目標 (包含 MyAlgo 1-6)
a.out: main.cpp $(dependency) $(algorithm) $(routing)
	g++ -fopenmp $(gcc_flag) main.cpp $(dependency) $(algorithm) $(routing) $(gurobi_flag) -o a.out
	> ../data/log/ZFA2_Purification_Stats.txt
	> ../data/log/ZFA_routing_Purification_Stats.txt
	./a.out

# --- 各個物件檔案 (Object Files) 編譯規則 ---

# 網路核心組件
Network/Node/Node.o: Network/Node/Node.cpp Network/Node/Node.h config.h
	g++ -c $(gcc_flag) Network/Node/Node.cpp -o Network/Node/Node.o

Network/Shape/Shape.o: Network/Shape/Shape.cpp Network/Shape/Shape.h Network/Purification/Purification.h config.h
	g++ -c $(gcc_flag) Network/Shape/Shape.cpp -o Network/Shape/Shape.o

Network/Graph/Graph.o: Network/Graph/Graph.cpp Network/Graph/Graph.h Network/Purification/Purification.h config.h
	g++ -c $(gcc_flag) Network/Graph/Graph.cpp -o Network/Graph/Graph.o

Algorithm/AlgorithmBase/AlgorithmBase.o: Algorithm/AlgorithmBase/AlgorithmBase.h Algorithm/AlgorithmBase/AlgorithmBase.cpp config.h
	g++ -c $(gcc_flag) Algorithm/AlgorithmBase/AlgorithmBase.cpp -o Algorithm/AlgorithmBase/AlgorithmBase.o

# MyAlgo 1-6 編譯規則
Algorithm/MyAlgo1/MyAlgo1.o: Algorithm/MyAlgo1/MyAlgo1.cpp Algorithm/MyAlgo1/MyAlgo1.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo1/MyAlgo1.cpp -o Algorithm/MyAlgo1/MyAlgo1.o
Algorithm/MyAlgo2/MyAlgo2.o: Algorithm/MyAlgo2/MyAlgo2.cpp Algorithm/MyAlgo2/MyAlgo2.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo2/MyAlgo2.cpp -o Algorithm/MyAlgo2/MyAlgo2.o
Algorithm/MyAlgo3/MyAlgo3.o: Algorithm/MyAlgo3/MyAlgo3.cpp Algorithm/MyAlgo3/MyAlgo3.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo3/MyAlgo3.cpp -o Algorithm/MyAlgo3/MyAlgo3.o
Algorithm/MyAlgo4/MyAlgo4.o: Algorithm/MyAlgo4/MyAlgo4.cpp Algorithm/MyAlgo4/MyAlgo4.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo4/MyAlgo4.cpp -o Algorithm/MyAlgo4/MyAlgo4.o
Algorithm/MyAlgo5/MyAlgo5.o: Algorithm/MyAlgo5/MyAlgo5.cpp Algorithm/MyAlgo5/MyAlgo5.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo5/MyAlgo5.cpp -o Algorithm/MyAlgo5/MyAlgo5.o
Algorithm/MyAlgo6/MyAlgo6.o: Algorithm/MyAlgo6/MyAlgo6.cpp Algorithm/MyAlgo6/MyAlgo6.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo6/MyAlgo6.cpp -o Algorithm/MyAlgo6/MyAlgo6.o

# Werner 演算法編譯規則
Algorithm/WernerAlgo/WernerAlgo.o: Algorithm/WernerAlgo/WernerAlgo.cpp Algorithm/WernerAlgo/WernerAlgo.h config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo/WernerAlgo.cpp -o Algorithm/WernerAlgo/WernerAlgo.o
Algorithm/WernerAlgo2/WernerAlgo2.o: Algorithm/WernerAlgo2/WernerAlgo2.cpp Algorithm/WernerAlgo2/WernerAlgo2.h Network/Purification/Purification.h config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo2/WernerAlgo2.cpp -o Algorithm/WernerAlgo2/WernerAlgo2.o
Algorithm/WernerAlgo_routing/WernerAlgo_routing.o: Algorithm/WernerAlgo_routing/WernerAlgo_routing.cpp Algorithm/WernerAlgo_routing/WernerAlgo_routing.h Network/Purification/Purification.h config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo_routing/WernerAlgo_routing.cpp -o Algorithm/WernerAlgo_routing/WernerAlgo_routing.o
Algorithm/WernerAlgo3/WernerAlgo3.o: Algorithm/WernerAlgo3/WernerAlgo3.cpp Algorithm/WernerAlgo3/WernerAlgo3.h config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo3/WernerAlgo3.cpp -o Algorithm/WernerAlgo3/WernerAlgo3.o
Algorithm/WernerAlgo_UB/WernerAlgo_UB.o: Algorithm/WernerAlgo_UB/WernerAlgo_UB.cpp Algorithm/WernerAlgo_UB/WernerAlgo_UB.h config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo_UB/WernerAlgo_UB.cpp -o Algorithm/WernerAlgo_UB/WernerAlgo_UB.o
Algorithm/GreedySP/GreedySP.o: Algorithm/GreedySP/GreedySP.cpp Algorithm/GreedySP/GreedySP.h Network/Purification/Purification.h config.h
	g++ -c $(gcc_flag) Algorithm/GreedySP/GreedySP.cpp -o Algorithm/GreedySP/GreedySP.o

# 時間量測版本演算法編譯規則
Algorithm/WernerAlgo_time/WernerAlgo_time.o: Algorithm/WernerAlgo_time/WernerAlgo_time.cpp Algorithm/WernerAlgo_time/WernerAlgo_time.h config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo_time/WernerAlgo_time.cpp -o Algorithm/WernerAlgo_time/WernerAlgo_time.o
Algorithm/WernerAlgo2_time/WernerAlgo2_time.o: Algorithm/WernerAlgo2_time/WernerAlgo2_time.cpp Algorithm/WernerAlgo2_time/WernerAlgo2_time.h config.h
	g++ -c $(gcc_flag) Algorithm/WernerAlgo2_time/WernerAlgo2_time.cpp -o Algorithm/WernerAlgo2_time/WernerAlgo2_time.o
Algorithm/MyAlgo1_time/MyAlgo1_time.o: Algorithm/MyAlgo1_time/MyAlgo1_time.cpp Algorithm/MyAlgo1_time/MyAlgo1_time.h config.h
	g++ -c $(gcc_flag) Algorithm/MyAlgo1_time/MyAlgo1_time.cpp -o Algorithm/MyAlgo1_time/MyAlgo1_time.o

# 路由與路徑方法編譯規則 (REPS 需連結 Gurobi)
Network/PathMethod/REPS/REPS.o: Network/PathMethod/REPS/REPS.cpp Network/PathMethod/REPS/REPS.h config.h
	g++ -c $(gcc_flag) Network/PathMethod/REPS/REPS.cpp $(gurobi_flag) -o Network/PathMethod/REPS/REPS.o
Network/PathMethod/Greedy/Greedy.o: Network/PathMethod/Greedy/Greedy.cpp Network/PathMethod/Greedy/Greedy.h config.h
	g++ -c $(gcc_flag) Network/PathMethod/Greedy/Greedy.cpp -o Network/PathMethod/Greedy/Greedy.o
Network/PathMethod/QCAST/QCAST.o: Network/PathMethod/QCAST/QCAST.cpp Network/PathMethod/QCAST/QCAST.h config.h
	g++ -c $(gcc_flag) Network/PathMethod/QCAST/QCAST.cpp -o Network/PathMethod/QCAST/QCAST.o
Network/PathMethod/PathMethodBase/PathMethod.o: Network/PathMethod/PathMethodBase/PathMethod.cpp Network/PathMethod/PathMethodBase/PathMethod.h config.h
	g++ -c $(gcc_flag) Network/PathMethod/PathMethodBase/PathMethod.cpp -o Network/PathMethod/PathMethodBase/PathMethod.o

config.o: config.h config.cpp
	g++ -c $(gcc_flag) config.cpp -o config.o

# --- 清除工具 ---
clean:
	rm -f a.out zfa_run zfa2_run runtime
	find . -type f -name '*.o' -delete
	find . -type f -name '*.h.gch' -delete
	rm -f ../data/ans/*.ans
	rm -f ../data/log/*
	rm -f ../data/input/*

.PHONY: clean all zfa_run zfa2_run runtime
