Online Spatiotemporal Routing for Congestion Minimization
========================================================================

This repository contains the source codes and data for this paper. 

Usage
---------------

### Environment

g++ version: 7.4.0 

OS: Ubuntu/CentOS (Unix)

### Compilation

cd code

make

### Execution

In the folder code/,

./shortesPath [network path] [query path]

./IPM [network path] [query path]

./BR [network path] [query path] [detour factor]

./GRO [network path] [query path] [detour factor]

./SOR [network path] [query path] [detour factor]

./SRH [network path] [query path] [detour factor] [candidate path]

For dynamic edge weights,

./shortesPathD [network path] [query path] [updatedEdges path]

./SORD [network path] [query path] [detour factor] [updatedEdges path]

./SRHD [network path] [query path] [detour factor] [candidate path] [updatedEdges path]

The default setting uses

network path=../data/network_NYC.txt

query_path=../data/qnNYC/qn10000

detour factor=0.05

candidate path=../data/candidate5_NYC.txt

updatedEdges path=../data/update/update1.txt

A complete statement is ./SRH ../data/network_NYC.txt ../data/qnNYC/qn5000 0.05 ../data/candidate5_NYC.txt

The output of the screen contains the maximum value, execution time, and memory cost and a detailed load distribution recorded in load_[SHORT/IPM/BRSOR/SRH].txt, where each line contains the edge id, the time step, and the corresponding load. If the dynamic versions of algorithms are executed, there is an update time following the execution time.

You may want to generate the candidate yourself by Statistic.

./Statistic [network path] [history path] [starting step] [ending step] [output filename] [detour factor]

A default run is ./Statistic ../data/network_NYC.txt ../history/NYC/ 1020 1080 candidate.txt 0.05

Note that 1020 stands for 17:00 since one step corresponds to one minute.

### Data Description

The data used in the experiments are stored in data/qnNYC and data/qnBJ. 

Each network file contains the number of nodes n and edges m in the first line. The following m lines have the origin, the destination node, and the weights of edges w_e. Note that one unit of edge weights is one second.

Each query file contains the number of queries |Q|, the number of steps T, and the detour factor a in the first line, and the step t, the origin node, the destination node of each query in the following |Q| lines.

Each candidate file contains the number of candidate edge-step pair |C| and the starting step in the first line, the mean value and the confidence radius of each Q^t in the following T lines, and the edge id and the step id (e,\tau) in the following |C|.

Each history file in the history folder contains the queries in the previous 20 days in the same format as a query file.

Each update file in the update folder contains the number of edge weight updates in the first line, the edge id, the new weights, and the update time in the following lines.


