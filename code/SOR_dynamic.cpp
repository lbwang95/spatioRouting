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
typedef pair<int, double> ID;
typedef pair<double, ID> Tri;
const int MAX_V = 270000, MAX_M = 740000, MAX_T = 1540, MAX_Q = 30000;// max number of nodes, edges, and steps
const double EPS = 1e-8;
int pre[MAX_V], pre_i[MAX_V], V, M, T;
int nquery, ndnet, nslots;
double cap[MAX_M];
double load[MAX_M][MAX_T], max_load;
double detour, alpha = 0.01, sdist[MAX_V], cong[MAX_V], Lambda = 1;
DDV shortest, c_short;
set<Tri> trav[MAX_M];//for the index trav(e)
struct edge{
    int to, id;
    double cost;
    edge(int _to, int _id, double _cost){
        to = _to;
        id = _id;
        cost = _cost;
    }
};

vector<edge> G1[MAX_V], G2[MAX_V];
vector<int> pathset[MAX_Q];
int pathnum;
int mod = 30000;

//calculate the load after the path is recommended
double cal_load(int t, int s, int d, vector<int> &path){
    int v = s;
    int j = 0, tau = t;
    double ub = 0, theta = 0;
    pathset[pathnum%mod].push_back(v);
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        pathset[pathnum%mod].push_back(path[i]);
        int flag = 0;
        for (; j < ub; j++){
            load[e.id][tau + j] += 1 / cap[e.id];
            double x_variable = pow((1 + 1 / (2 * Lambda)), load[e.id][tau + j]) *alpha/ M;
            theta += x_variable;
            if(x_variable>exp(0.5))
                Lambda *= 2;
            max_load = max(max_load, load[e.id][tau + j]);
            if (flag == 0){
                flag = 1;
                trav[e.id].insert(Tri(t+ub, ID(pathnum%mod, t)));
            }
        }
        v = e.to;
    }
    pathnum++;
    if(theta>Lambda)
        Lambda *= 2;
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
//compute the metric value for a path
double update(vector<int> &path, int t, int s, int d){
    int v = s;
    double ub = 0, conE = 0;
    int p = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; p < ub; p++){
            conE += pow((1 + 1 / (2 * Lambda)), load[e.id][t + p]);
        }
        v = e.to;
    }
    return conE;
}
//first do a backward search on the reverse graph G2
double dijkstra2(int t, int s, int d){
    priority_queue<P, vector<P>, greater<P> > que;
    for (int i = 0; i <= V;i++)
        sdist[i] = DBL_MAX;
    sdist[d] = 0;
    que.push(P(0, d));
    while (!que.empty()){
        P p = que.top();
        que.pop();
        int v = p.second;
        if (sdist[v] < p.first)
            continue;
        if(v==s){
            vector<int> path;
            for (; v!= d; v=pre[v]){
                int eid = G2[pre[v]][pre_i[v]].id;
                for (int i = 0; i < G1[v].size();i++){
                    if(G1[v][i].id==eid)
                        path.push_back(i);
                }
            }
            shortest.second.first = sdist[s];
            shortest.second.second = path;
            shortest.first = update(path, t, s, d);
            return sdist[s];
        }
        for (int i = 0; i < G2[v].size(); i++) {
            edge e = G2[v][i];
            if (sdist[e.to] > sdist[v] + e.cost){
                sdist[e.to] = sdist[v] + e.cost;
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(P(sdist[e.to], e.to));
            }
        }
    }       
    return sdist[s];
}
//find the path with the minimum metric value
void c_dijkstra(int t, int s, int d){
    fill(pre, pre + V, -1);//pre store the previous node in the path
    fill(pre_i, pre_i + V, -1);//pre_i store the index of the previous node in the path
    priority_queue<DDI, vector<DDI>, comp > que;
    for (int i = 0; i <= V;i++){
        cong[i] = DBL_MAX;
    }    
    cong[s] = 0;
    que.push(DDI(0, P(0, s)));
    while (!que.empty()){
        DDI p = que.top();
        que.pop();
        double congv = p.first;
        int v = p.second.second;
        double dis = p.second.first;
        if(v==d){
            vector<int> path;
            for (int v=d; v!= s; v=pre[v])
                path.push_back(pre_i[v]);
            reverse(path.begin(), path.end());
            c_short.second.first = dis;
            c_short.second.second = path;
            c_short.first = update(path, t, s, d);
        }
        if (cong[v] < congv)
            continue;
        for (int i = 0; i < G1[v].size(); i++) {
            edge e = G1[v][i];
            if (dis + e.cost + sdist[e.to] > (1 + detour) * sdist[s])
                continue;
            double ub =  dis + e.cost, conE = 0;
            for (int j = ceil(dis); j < ub; j++){
                conE += pow((1 + 1 / (2 * Lambda)), load[e.id][t + j]);
            }
            if(cong[e.to] > congv + conE){
                cong[e.to] = congv + conE;
                pre[e.to] = v;
                pre_i[e.to] = i;
                que.push(DDI(cong[e.to], P(ub, e.to)));
            }
        }
    }
}

void SOR(int t, int s, int d){
    if (sdist[s] == DBL_MAX){
        printf("!errrrrrrrrrrrrrrrrrrrrrrrrr:%d %d %d!\n", t, s, d);
        return;
    }
    double shortestdist = shortest.second.first;
    if(shortest.first<cong[d]+EPS){
        //printf("&\n");
        cal_load(t, s, d, shortest.second.second);
        #ifdef WATCH_MEM
            watchSolutionOnce(getpid(), usedMemory);
        #endif
        return;
    }
    if(c_short.second.first<(1 + detour) * shortestdist+EPS){
        //printf("#\n");
        cal_load(t, s, d, c_short.second.second);
        #ifdef WATCH_MEM
            watchSolutionOnce(getpid(), usedMemory);
        #endif        
        return;
    }
    printf("%d %d %d:%f %f %f %f ", t, s, d, cong[d], shortest.first, c_short.second.first, sdist[s]);
    return;
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
                neworigin = e.to;
            }
            int flag = 0;
            for (; j < ub; j++){
                load[e.id][tau + j] -= 1 / cap[e.id];
                if (flag == 0){
                    flag = 1;
                    trav[e.id].erase(Tri(tau+ub, ID(pathid, deptime)));
                }
            }
            v = e.to;
        }
        dijkstra2(updatetime, neworigin, v);
        c_dijkstra(updatetime, neworigin, v);
        SOR(updatetime, neworigin, v);
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
    fscanf(fp_network, "%d%d", &V, &M);
    int a, b, t;
    double c;
    for (int i = 0; i < M; i++) {
        fscanf(fp_network,"%d%d%lf", &a, &b, &c);
        alledges.push_back(PI(a, b));
        G1[a].push_back(edge(b, i, c / 60));
        G2[b].push_back(edge(a, i, c / 60));
        cap[i] = c*40000/3600/50*2;
        if(cap[i]<2)
            cap[i] = 100;
    }
    fscanf(fp_query, "%d%d%lf", &nquery, &T, &detour);
    if (argc > 3)
        detour = stod(string(argv[3]));
    else
        detour = 0.05;
    if (argc > 4)
        fp_update = fopen(argv[4], "r");
    else
        fp_update=fopen("../data/update/update1.txt","r");
    int nupdates;
    fscanf(fp_update,"%d", &nupdates);
    for (int i = 0; i < nupdates;i++){
        fscanf(fp_update, "%d%lf%d", &a, &c, &t);
        updates[t].push_back(P(c, a));
    }
    clock_t proc_s = clock();
    double updateTimeCost = 0;
    int _t = -1;
    for (int i = 0; i < nquery;i++){
        fscanf(fp_query, "%d%d%d", &t, &a, &b);
        if(a==b)
            continue;
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
        dijkstra2(t, a, b);
        c_dijkstra(t, a, b);
        SOR(t, a, b);
    }
    clock_t proc_e = clock();
    double proc_d = proc_e - proc_s;
    #ifdef WATCH_MEM
        printf("%f %f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, (double)updateTimeCost / CLOCKS_PER_SEC, usedMemory / 1024.0);
    #else
        printf("%f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, (double)updateTimeCost / CLOCKS_PER_SEC);
    #endif

    FILE *fp = fopen("load_SORD.txt", "w");
    for (int i = 0; i < M;i++)
        for (int j = 0; j < T+40;j++){
            if (load[i][j] > max_load/3*2){
                fprintf(fp, "%d %d %f\n", i, j, load[i][j]);
            }
        }
    fclose(fp);
    return 0;
}
