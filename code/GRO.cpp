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
double load[MAX_M][MAX_T], vol[MAX_M][MAX_T], max_load, dist[MAX_V], loade[MAX_M];
double detour, alpha = 0.01, Lambda = 1;
DDV shortest, c_short;
double N;
vector<PI> queryset[1000];

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
    int v = s, j = 0, tau = t;
    double ub = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; j < ub; j++){
            load[e.id][tau + j] += 1 / cap[e.id];
            vol[e.id][tau + j] += 1;
            max_load = max(max_load, load[e.id][tau + j]);
        }
        v = e.to;
    }
    return max_load;
}

double WT(int t, int s, int d, vector<int> &path){
    int v = s, j = 0, tau = t;
    double ub = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; j < ub; j++){
            vol[e.id][tau + j] += 1;
        }
        v = e.to;
    }
}

double traveltimel(double curt,int eid, double oritt){
    int curet = ceil(curt);
    return oritt * (1 + 0.15 * pow((double)load[eid][curet], 4));
}

vector<vector<int>> pathsetl, pathsetv;
void dijkstral(int t, int s, int d){
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
            double tt = traveltimel(dist[v], e.id, e.cost);
            if (dist[e.to] > dist[v] + tt){
                dist[e.to] = dist[v] + tt;
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    if(dist[d]==DBL_MAX)
        printf("errrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr\n");
    WT(t, s, d, path);
    pathsetl.push_back(path);
    #ifdef WATCH_MEM
        watchSolutionOnce(getpid(), usedMemory);
    #endif
}

double traveltimev(double curt,int eid, double oritt){
    int curet = ceil(curt);
    return oritt * (1 + 0.15 * pow(vol[eid][curet]/cap[eid], 4));
}

void dijkstrav(int t, int s, int d){
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
            double tt = traveltimev(dist[v], e.id, e.cost);
            if (dist[e.to] > dist[v] + tt){
                dist[e.to] = dist[v] + tt;
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(P(dist[e.to], e.to));
            }
        }
    }
    if(dist[d]==DBL_MAX)
        printf("errrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr\n");
    pathsetv.push_back(path);
    #ifdef WATCH_MEM
        watchSolutionOnce(getpid(), usedMemory);
    #endif
}

double RWT(int t, int s, int d, vector<int> &path){
    int v = s, j = 0, tau = t;
    double ub = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; j < ub; j++){
            vol[e.id][tau + j] -= 1;
        }
        v = e.to;
    }
    return max_load;
}

int prev_t;
void GRO(int t, int s, int d){
    if (t != prev_t){
        prev_t = t;
        for (int i = 0; i < queryset[t-1].size();i++){
            int ori = queryset[t - 1][i].first, des = queryset[t - 1][i].second;
            dijkstrav(t - 1, ori, des);
        }
        int nl = pathsetl.size();
        for (int i = nl - 1; i >= 0; i--){
            int ori = queryset[t - 1][i].first, des = queryset[t - 1][i].second;
            RWT(t - 1, ori, des, pathsetl[i]);
            pathsetl.pop_back();
        }
        int nv = pathsetv.size();
        for (int i = nv - 1; i >= 0; i--){
            int ori = queryset[t - 1][i].first, des = queryset[t - 1][i].second;
            cal_load(t - 1, ori, des, pathsetv[i]);
            pathsetv.pop_back();
        }
    }
    dijkstral(t, s, d);
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
        queryset[t].push_back(PI(a, b));
        N = i + 1;
        //if(a==b)
        //    continue;
        GRO(t, a, b);
    }
    GRO(prev_t + 1, 0, 0);
    clock_t proc_e = clock();
    double proc_d = proc_e - proc_s;
    #ifdef WATCH_MEM
        printf("%f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, usedMemory / 1024.0);
    #else
        printf("%f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC);
    #endif

    FILE *fp = fopen("load_GRO.txt", "w");
    for (int i = 0; i < M;i++)
        for (int j = 0; j < T+40;j++){
            if (load[i][j] > max_load/3*2){
                fprintf(fp, "%d %d %f\n", i, j, load[i][j]);
            }
        }
    return 0;
}
