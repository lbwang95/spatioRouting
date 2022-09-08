#include<bits/stdc++.h>
#ifdef WATCH_MEM
	#include "monitor.h"
    int usedMemory;
#endif
using namespace std;
typedef pair<double, int> P;
typedef pair<int, int> PI;
typedef pair<int, double> ID;
typedef pair<double, ID> Tri;
const int MAX_V = 270000, MAX_M = 740000, MAX_T = 1540, MAX_Q = 30000;//max number of nodes, edges, and steps
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
int pre[MAX_V], pre_i[MAX_V], V, M, T, nquery;
double load[MAX_M][MAX_T], max_load;
double detour, dist[MAX_V];
double cap[MAX_M];
set<Tri> trav[MAX_M];//for the index trav(e)
vector<int> pathset[MAX_Q];
int pathnum;
int mod = 30000;
//calculate the load after the path is recommended
double cal_load(int t, int s, int d, vector<int> path){
    int v = s, j = 0, tau = t;
    double ub = 0;
    pathset[pathnum%mod].push_back(v);
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        pathset[pathnum%mod].push_back(path[i]);
        int flag = 0;
        for (; j < ub; j++){
            load[e.id][tau + j] += 1 / cap[e.id];
            max_load = max(max_load, load[e.id][tau + j]);
            if (flag == 0){
                flag = 1;
                trav[e.id].insert(Tri(t+ub, ID(pathnum%mod, t)));
            }
        }
        v = e.to;
    }
    pathnum++;
    return max_load;
}

vector<int> dijkstra(int t, int s, int d){
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
    //printf("%d\n", cal_load(t, s, d, path));
    cal_load(t, s, d, path);
    #ifdef WATCH_MEM
        watchSolutionOnce(getpid(), usedMemory);
    #endif
    return path;
}

void edgeWeightUpdate(int s, int d, double newcost, double updatetime){
    int ti = -1;
    for (int i = 0; i < G1[s].size();i++){
        if (G1[s][i].to == d){
            ti = i;
            break;
        }
    }
    if (ti == -1)
        return;
    int eid = G1[s][ti].id;
    double oldcost = G1[s][ti].cost;
    G1[s][ti].cost = newcost;
    if (abs(newcost-oldcost)<1e-5)
        return;
    set<Tri>:: iterator iter;
    iter = trav[eid].lower_bound(Tri(updatetime, ID(-1, 0)));
    trav[eid].erase(trav[eid].begin(),iter);
    vector<Tri> travcopy;
    for(set<Tri>::iterator it1=trav[eid].begin(); it1!=trav[eid].end(); ++it1)
        travcopy.push_back(*it1);
    for (int i = 0; i < travcopy.size();i++){
        Tri tmp = travcopy[i];
        int pathid=tmp.second.first;
        double arrtime = tmp.first;
        int deptime = tmp.second.second;
        int v = pathset[pathid][0], j = 0, tau = deptime, neworigin=-1;
        double ub = 0, theta = 0, newDeptime;
        for(int i = 1; i < pathset[pathid].size(); i++){
            edge e = G1[v][pathset[pathid][i]];
            if(e.id==eid)
                ub += oldcost;
            else
                ub += e.cost;
            if (tau + ub > updatetime && neworigin == -1){
                newDeptime = tau + ub + (tau + ub - updatetime) * newcost / oldcost;
                neworigin = v;
            }
            int flag = 0;
            for (; j < ub; j++){
                if (j >= updatetime)
                    load[e.id][tau + j] -= 1 / cap[e.id];
                if (flag == 0){
                    flag = 1;
                    trav[e.id].erase(Tri(tau+ub, ID(pathid, deptime)));
                }
            }
            v = e.to;
        }
        dijkstra(updatetime, neworigin, v);
    }
    
}
vector<P> updates[1540];
vector<PI> alledges;
int main(int argc , char * argv[]){
    FILE *fp_query, *fp_network, *fp_update;
    if (argc > 1)
        fp_network = fopen(argv[1], "r");
    else
        fp_network = fopen("../data/network_NYC.txt", "r");
    if (argc > 2)
        fp_query = fopen(argv[2], "r");
    else
        fp_query = fopen("../data/syn/qn/qn20000", "r");
    fscanf(fp_network,"%d%d", &V, &M);
    int a, b, t;
    double c;
    for (int i = 0; i < M; i++) {
        fscanf(fp_network,"%d%d%lf", &a, &b, &c);
        alledges.push_back(PI(a, b));
        G1[a].push_back(edge(b, i, c/60));
        cap[i] = c*40000/3600/50*2;
        if(cap[i]<2)
            cap[i] = 100;
    }
    if (argc > 3)
        fp_update = fopen(argv[3], "r");
    else
        fp_update=fopen("../data/update/update1.txt","r");
    int nupdates;
    fscanf(fp_update,"%d", &nupdates);
    for (int i = 0; i < nupdates;i++){
        fscanf(fp_update, "%d%lf%d", &a, &c, &t);
        updates[t].push_back(P(c, a));
    }
    double updateTimeCost = 0;
    int _t = -1;
    clock_t proc_s = clock();
    fscanf(fp_query, "%d%d%lf", &nquery, &T, &detour);
    for (int i = 0; i < nquery;i++){
        fscanf(fp_query, "%d%d%d", &t, &a, &b);
        if (t != _t){
            _t = t;
            clock_t proc_sup = clock();
            for (int j = 0; j < updates[t].size();j++){
                PI eab = alledges[updates[t][j].second];
                edgeWeightUpdate(eab.first, eab.second, updates[t][j].first, t);
            }
            clock_t proc_eup = clock();
            updateTimeCost += proc_eup - proc_sup;
        }
        dijkstra(t, a, b);
    }
    clock_t proc_e = clock();
    double proc_d = proc_e - proc_s;
    #ifdef WATCH_MEM
        printf("%f %f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, (double)updateTimeCost / CLOCKS_PER_SEC, usedMemory / 1024.0);
    #else
        printf("%f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, (double)updateTimeCost / CLOCKS_PER_SEC);
    #endif

    FILE *fp = fopen("load_SHORTD.txt", "w");
    for (int i = 0; i < M;i++)
        for (int j = 0; j < T+40;j++){
            if (load[i][j] > 0)//max_load/3*2)            
	    {
                fprintf(fp, "%d %d %f\n", i, j, load[i][j]);
            }
        }
    fclose(fp);
    return 0;
}
