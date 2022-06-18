#include<bits/stdc++.h>
#ifdef WATCH_MEM
	#include "monitor.h"
    int usedMemory;
#endif
using namespace std;
typedef pair<double, int> P;
const int MAX_V = 270000, MAX_M = 740000, MAX_T = 1540;//max number of nodes, edges, and steps
struct edge{
    int to, id;
    double cost, ocost;
    edge(int _to, int _id, double _cost){
        to = _to;
        id = _id;
        cost = _cost;
        ocost = _cost;
    }
};

vector<edge> G1[MAX_V];
int pre[MAX_V], pre_i[MAX_V], V, M, T, nquery;
double load[MAX_M][MAX_T], max_load;
double detour, dist[MAX_V];
double cap[MAX_M];

//calculate the load after the path is recommended
double cal_load(int t, int s, int d, vector<int> path){
    int v = s, j = 0, tau = t;
    double ub = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.ocost;
        for (; j < ub; j++){
            load[e.id][tau + j] += 1 / cap[e.id];
            max_load = max(max_load, load[e.id][tau + j]);
        }
        v = e.to;
    }
    return max_load;
}
//compute the shortest path for the alternative path
vector<int> dijkstra1(int t, int s, int d){
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i <= V;i++)
        dist[i] = DBL_MAX;
    memset(pre, -1, sizeof(pre));//pre store the previous node in the path
    memset(pre_i, -1, sizeof(pre_i));//pre_i store the index of the previous node in the path
    dist[s] = 0;
    que.push(P(0, s));
    vector<int> path;
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (dist[v] < p.first)
            continue;
        if (v == d){
            for (; v != s; v = pre[v])
                path.push_back(pre_i[v]);
            reverse(path.begin(), path.end());
            break;
        }
        for (int i = 0; i < G1[v].size(); i++) {
            edge e = G1[v][i];
            if (dist[e.to] > dist[v] + e.cost){
                dist[e.to] = dist[v] + e.cost;
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    if(dist[d]==DBL_MAX)
        printf("errrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr\n");
    return path;
}
//compute the alternative path after the shortest path is computed
vector<int> dijkstra(int t, int s, int d){
    vector<int> path1 = dijkstra1(t, s, d);
    int u = s;
    for(int i = 0; i < path1.size(); i++){
        edge e = G1[u][path1[i]];
        if(cap[e.id]>=2)
            G1[u][path1[i]].cost *= 1.4;
        u = e.to;
    }
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i <= V;i++)
        dist[i] = DBL_MAX;
    memset(pre, -1, sizeof(pre));
    memset(pre_i, -1, sizeof(pre_i));
    dist[s] = 0;
    que.push(P(0, s));
    vector<int> path;
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (dist[v] < p.first)
            continue;
        if (v == d){
            for (; v != s; v = pre[v])
                path.push_back(pre_i[v]);
            reverse(path.begin(), path.end());
            break;
        }
        for (int i = 0; i < G1[v].size(); i++) {
            edge e = G1[v][i];
            if (dist[e.to] > dist[v] + e.cost){
                dist[e.to] = dist[v] + e.cost;
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    if(dist[d]==DBL_MAX)
        printf("errrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr\n");
    //printf("%d\n", cal_load(t, s, d, path));
    cal_load(t, s, d, path);
    u = s;
    for(int i = 0; i < path1.size(); i++){
        edge e = G1[u][path1[i]];
        if(cap[e.id]>=2)
            G1[u][path1[i]].cost = G1[u][path1[i]].ocost;
        u = e.to;
    }
    #ifdef WATCH_MEM
        watchSolutionOnce(getpid(), usedMemory);
    #endif
    return path;
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
    if (argc > 3)
        detour = stod(string(argv[3]));
    else
        detour = 0.05;
    fscanf(fp_network,"%d%d", &V, &M);
    int a, b, t;
    double c;
    for (int i = 0; i < M; i++) {
        fscanf(fp_network,"%d%d%lf", &a, &b, &c);
        G1[a].push_back(edge(b, i, c/60));
        cap[i] = c*40000/3600/50*2;
        if(cap[i]<2)
            cap[i] = 100;
    }
    clock_t proc_s = clock();
    fscanf(fp_query, "%d%d%lf", &nquery, &T, &detour);
    for (int i = 0; i < nquery;i++){
        fscanf(fp_query, "%d%d%d", &t, &a, &b);
        //printf("%d %d %d\n", t, a, b);
        dijkstra(t, a, b);
    }
    clock_t proc_e = clock();
    double proc_d = proc_e - proc_s;
    #ifdef WATCH_MEM
        printf("%f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, usedMemory / 1024.0);
    #else
        printf("%f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC);
    #endif

    FILE *fp = fopen("load_IPM.txt", "w");
    for (int i = 0; i < M;i++)
        for (int j = 0; j < T+40;j++){
            if (load[i][j] > max_load/3*2){
                fprintf(fp, "%d %d %f\n", i, j, load[i][j]);
            }
        }
    return 0;
}
