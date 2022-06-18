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
typedef pair<double, DV > DDV;
typedef pair<double, P> DDI;
const int MAX_V = 270000, MAX_M = 740000, MAX_T = 1540;//max number of nodes, edges, and steps
const double EPS = 1e-8;
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
double Q[MAX_T][2];
set<PI> cands;
int pre[MAX_V], pre_i[MAX_V], V, M, T, C;
int nquery, ndnet, nslots;
double cand[MAX_M][MAX_T];
double cap[MAX_M];
double load[MAX_M][MAX_T], max_load;
double detour, alpha = 0.01, sdist[MAX_V], cong[MAX_V], Lambda = 1;
DDV shortest, c_short;
double bal = 0.001;

//calculate the load after the path is recommended
double cal_load(int t, int s, int d, vector<int> path){
    int v = s, j = 0, tau = t;
    double ub = 0, theta = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; j < ub; j++){
            load[e.id][tau + j] += 1 / cap[e.id];
            double x_variable = (cand[e.id][tau + j]+bal)*pow((1 + 1 / (2 * Lambda)), load[e.id][tau + j]) *alpha/ M;
            theta += x_variable;
            if(x_variable>exp(0.5))
                Lambda *= 2;
            max_load = max(max_load, load[e.id][tau + j]);
        }
        v = e.to;
    }
    if(theta>Lambda)
        Lambda *= 2;
    return max_load;
}
//compute the metric value for a path
double update(vector<int> &path, int t, int s, int d){
    double ub = 0, conE = 0;
    int v = s, p = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; p < ub; p++){
            conE += (cand[e.id][t+p]+bal)*pow((1 + 1 / (2 * Lambda)), load[e.id][t + p]);
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
                pre[e.to] = v;
                pre_i[e.to] = i;
                sdist[e.to] = sdist[v] + e.cost;
                que.push(P(sdist[e.to], e.to));
            }
        }
    }       
    return sdist[s];
}

struct comp {
    bool operator()(DDI& a, DDI& b) {
        if (a.first == b.first)
            return a.second.first > b.second.first;
        else
            return a.first > b.first;
    }
};
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
                conE += (cand[e.id][t+j]+bal) * pow((1 + 1 / (2 * Lambda)), load[e.id][t + j]);
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


void SRH(int t, int s, int d){
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

int main(int argc , char * argv[]){
    FILE *fp_query, *fp_network, *fp_candidate;
    if (argc > 1)
        fp_network = fopen(argv[1], "r");
    else
        fp_network = fopen("../data/network_NYC.txt", "r");
    if (argc > 2)
        fp_query = fopen(argv[2], "r");
    else
        fp_query = fopen("../data/syn/qn/qn20000", "r");
    if (argc > 4)
        fp_candidate = fopen(argv[4], "r");
    else
        fp_candidate=fopen("../data/syn/candidate.txt","r");
    fscanf(fp_network, "%d%d", &V, &M);
    int a, b, t;
    double c;
    for (int i = 0; i < M; i++) {
        fscanf(fp_network,"%d%d%lf", &a, &b, &c);
        G1[a].push_back(edge(b, i, c/60));
        G2[b].push_back(edge(a, i, c/60));
        cap[i] = c*40000/3600/50*2;
        if(cap[i]<2)
            cap[i] = 100;
    }
    fscanf(fp_query, "%d%d%lf", &nquery, &T, &detour);
    if (argc > 3)
        detour = stod(string(argv[3]));
    else
        detour = 0.05;
    int QT, offset;
    fscanf(fp_candidate,"%d%d%d", &C, &offset, &QT);
    for(int i = 0; i < QT; i++){
        fscanf(fp_candidate,"%lf%lf", &Q[i][0], &Q[i][1]);
    }
    for(int i = 0; i < C; i++){
        fscanf(fp_candidate,"%d%d%lf", &a, &b, &c);
        b -= offset;
        cand[a][b] = c;//1;
        cands.insert(PI(a, b));
    }
    clock_t proc_s = clock();
    int _t = -1;
    for (int i = 0; i < nquery;i++){
        fscanf(fp_query, "%d%d%d", &t, &a, &b);
        if(a==b)
            continue;
        vector<PI> tmpcand;
        if (t != _t){
            set<PI>::iterator itset;
            for (itset = cands.begin(); itset != cands.end();itset++){
                int e = itset->first, tau = itset->second;
                if(load[e][t]+Q[t][0]+Q[t][1]<max_load){
                    tmpcand.push_back(PI(e, t));
                    printf("|#####%f###%d####%f####%f####|", load[e][t], t,Q[t][0], max_load);
                }
            }
            for (int j = 0; j < tmpcand.size();j++)
                cand[tmpcand[j].first][tmpcand[j].second] = 0;
            _t = t;
        }
        dijkstra2(t, a, b);
        c_dijkstra(t, a, b);
        SRH(t, a, b);
        for (int j = 0; j < tmpcand.size();j++)
            cand[tmpcand[j].first][tmpcand[j].second] = 1;
    }
    clock_t proc_e = clock();
    double proc_d = proc_e - proc_s;
    #ifdef WATCH_MEM
        printf("%f %f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC, usedMemory / 1024.0);
    #else
        printf("%f %f\n", max_load, (double)proc_d / CLOCKS_PER_SEC);
    #endif

    FILE *fp = fopen("load_SRH.txt", "w");
    for (int i = 0; i < M;i++)
        for (int j = 0; j < T+40;j++){
            if(load[i][j]>max_load/3*2){
                fprintf(fp, "%d %d %f\n", i, j, load[i][j]);
            }
        }
    return 0;
}
