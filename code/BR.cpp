#include<bits/stdc++.h>
#ifdef WATCH_MEM
	#include "monitor.h"
    int usedMemory;
#endif
using namespace std;
typedef pair<int, int> PI;
typedef pair<double, int> P;
typedef pair<double, PI> DII;
typedef pair<double, vector<int>> DV;
typedef pair<int, vector<int>> IV;
typedef pair<double, IV> DIV;
typedef pair<double, DV > DDV;
typedef pair<double, P> DDI;
const int MAX_V = 270000, MAX_M = 740000, MAX_T = 1540;//max number of nodes, edges, and steps
const double EPS = 1e-8;
int pre[MAX_V], pre_i[MAX_V], V, M, T;
int nquery, ndnet, nslots;
double cap[MAX_M];
double load[MAX_M][MAX_T], max_load, dist[MAX_V], loade[MAX_M];
double detour, alpha = 0.01, Lambda = 1;
DDV shortest, c_short;
double N;

struct edge{
    int to, id;
    double cost;
    edge(int _to, int _id, double _cost){
        to = _to;
        id = _id;
        cost = _cost;
    }
};

vector<edge> G1[MAX_V];

//calculate the load after the path is recommended
double cal_load(int t, int s, int d, vector<int> &path){
    int v = s;
    int j = 0, tau = t;
    double ub = 0, theta = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; j < ub; j++){
            load[e.id][tau + j] += 1 / cap[e.id];
            max_load = max(max_load, load[e.id][tau + j]);
        }
        loade[e.id] += 1 / cap[e.id];
        v = e.to;
    }
    return max_load;
}

struct comp {
    bool operator()(DDI& a, DDI& b) {
        if (a.first == b.first)
            return a.second.first > b.second.first;
        else
            return a.first > b.first;
    }
};
vector<PI> resE;
//compute the metric value of a path based on the entropy
pair<double,double> update(vector<int> &path, int t, int s, int d){
    int v = s;
    double ub = 0, conE = 1;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        conE *= exp(-loade[e.id] / N * log((1 + loade[e.id]) / N));
        v = e.to;
        ub += e.cost;
    }
    return make_pair(ub, conE);
}
int flage[MAX_M];
//find the path which gives the minimum metric value
void BR(int t, int s, int d){
    double shortestdist, curconE = 10000000;
    vector<int> curpath;
    double penalty = 1.1;
    memset(flage, 0, sizeof(flage));
    for (int iter = 0; iter < 10; iter++){
        priority_queue<P, vector<P>, greater<P> > que;
        for (int i = 0; i <= V;i++)
            dist[i] = DBL_MAX;
        memset(pre, -1, sizeof(pre));//pre store the previous node in the path
        memset(pre_i, -1, sizeof(pre_i));//pre_i store the index of the previous node in the path
        dist[s] = 0;
        que.push(P(0, s));
        while (!que.empty()){
            P p = que.top();
            que.pop();
            int v = p.second;
            if (dist[v] < p.first)
                continue;
            if (v == d){
                vector<int> path;
                for (; v != s; v = pre[v]){
                    path.push_back(pre_i[v]);
                    flage[G1[pre[v]][pre_i[v]].id]++;
                }
                reverse(path.begin(), path.end());
                pair<double, double> tmp = update(path, t, s, d);
                if (iter == 0)
                    shortestdist = tmp.first;
                if (tmp.first < (1 + detour) * shortestdist && tmp.second < curconE){
                    curconE = tmp.second;
                    curpath = path;
                }
                break;
            }
            for (int i = 0; i < G1[v].size(); i++) {
                edge e = G1[v][i];
                if (dist[e.to] > dist[v] + e.cost * pow(penalty, flage[e.id])){
                    dist[e.to] = dist[v] + e.cost * pow(penalty, flage[e.id]);
                    pre[e.to] = v;
                    pre_i[e.to] = i;
                    que.push(P(dist[e.to], e.to));
                }
            }
        }
    }
    cal_load(t, s, d, curpath);
    #ifdef WATCH_MEM
        watchSolutionOnce(getpid(), usedMemory);
    #endif
    //printf("%d %d %d:%f\n", t, s, d, max_load);
    return;
}

int main(int argc , char * argv[]){
    FILE *fp_query, *fp_network;
    if (argc > 1)
        fp_network = fopen(argv[1], "r");
    else
        fp_network = fopen("../data/network_NYC.txt", "r");
    if (argc > 2)
        fp_query = fopen(argv[2], "r");
    else
        fp_query = fopen("../data/syn/qn/qn20000", "r");
    fscanf(fp_network, "%d%d", &V, &M);
    int a, b, t;
    double c;
    for (int i = 0; i < M; i++) {
        fscanf(fp_network,"%d%d%lf", &a, &b, &c);
        G1[a].push_back(edge(b, i, c / 60));
        cap[i] = c*40000/3600/50*2;
        if(cap[i]<2)
            cap[i] = 100;
    }
    clock_t proc_s = clock();
    fscanf(fp_query, "%d%d%lf", &nquery, &T, &detour);
    if (argc > 3)
        detour = stod(string(argv[3]));
    else
        detour = 0.05;
    for (int i = 0; i < nquery;i++){
        fscanf(fp_query, "%d%d%d", &t, &a, &b);
        N = i + 1;
        if(a==b)
            continue;
        BR(t, a, b);
    }
    clock_t proc_e = clock();
    double proc_d = proc_e - proc_s;
    #ifdef WATCH_MEM
        printf("%f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, usedMemory / 1024.0);
    #else
        printf("%f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC);
    #endif

    FILE *fp = fopen("load_BR.txt", "w");
    for (int i = 0; i < M;i++)
        for (int j = 0; j < T+40;j++){
            if (load[i][j] > max_load/3*2){
                fprintf(fp, "%d %d %f\n", i, j, load[i][j]);
            }
        }
    return 0;
}
