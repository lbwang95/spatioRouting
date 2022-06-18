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
const int MAX_V=270000,MAX_M=740000, MAX_T=1540, MAX_DAY=20;
int period = 5;
double det = 0.05;
double delta = 0.1, theta = 2, detour;
double cap[MAX_M];
map<int, map<int, double>> et2ha;
typedef struct candidate{
    int eid, t, nsamp;
    double rec[MAX_DAY];
    double ha, std, R, radius;
    candidate(int _eid, int _t, double _rec[MAX_DAY]){
        eid = _eid;
        t = _t;
        std = 0;
        R = 0;
        nsamp = MAX_DAY;
        double sum = 0;
        for (int i = 0; i < MAX_DAY;i++){
            rec[i] = _rec[i];
            sum += rec[i];
            R = max(R, (double)_rec[i]);
        }    
        ha = sum/MAX_DAY;
        et2ha[eid][t] = ha;
        for (int i = 0; i < MAX_DAY;i++)
            std += (rec[i] - ha) * (rec[i] - ha);
        std = sqrt(std/(MAX_DAY-1));
        radius = sqrt(R * R * log(2 / delta) / (2 * nsamp));
    }
}cand;

bool Cmp (cand a, cand b) {
    return a.ha + a.radius > b.ha + b.radius;
}

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
int V, M, T, C, pre[MAX_V], pre_i[MAX_V];
int nquery, ndnet, nslots;
map<int, map<int, double[MAX_DAY]> > sta;
int stat_Q[MAX_T][MAX_DAY];
double maxcong = 0;
vector<cand> sta_alg2;
map<int, int> e_alg3[MAX_T];
map<int, map<int, int>> t_alg3;
int t_lb = 1020, t_ub = t_lb + 60;
DDV shortest, c_short;
double alpha = 0.01, sdist[MAX_V], cong[MAX_V], Lambda = 1;
double load[MAX_M][MAX_T], max_load;
const double EPS = 1e-8;
//count the load after the path is recommended
double cal_load(int t, int s, int d, vector<int> &path){
    int v = s;
    int j = 0, tau = t;
    double ub = 0, theta = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; j < ub; j++){
            load[e.id][tau + j] += 1 / cap[e.id];
            double x_variable = pow((1 + 1 / (2 * Lambda)), load[e.id][tau + j]) *alpha/ M;
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
    fill(pre, pre + V, -1);
    fill(pre_i, pre_i + V, -1);
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
//count the number of trips 
void sor_stat(int day, int t, int s, int d){
    dijkstra2(t, s, d);
    c_dijkstra(t, s, d);
    vector<int> path;
    if (sdist[s] == DBL_MAX){
        printf("!errrrrrrrrrrrrrrrrrrrrrrrrr:%d %d %d!\n", t, s, d);
        return;
    }
    double shortestdist = shortest.second.first;
    if(shortest.first<cong[d]+EPS){
        //printf("&\n");
        path = shortest.second.second;
        cal_load(t, s, d, shortest.second.second);
    }
    else if(c_short.second.first<(1 + detour) * shortestdist+EPS){
        //printf("#\n");
        path=c_short.second.second;
        cal_load(t, s, d, c_short.second.second);
    }
    int v = s, j = 0, tau = t;
    double ub = 0;
    for(int i = 0; i < path.size(); i++){
        edge e = G1[v][path[i]];
        ub += e.cost;
        for (; j < ub; j++){
            if(sta.find(e.id)==sta.end()||sta[e.id].find((tau+j)/period)==sta[e.id].end())
            {
                for (int l = 0; l < MAX_DAY;l++)
                    sta[e.id][(tau + j)/period][l] = 0;
                sta[e.id][(tau + j) / period][day] = 1 / cap[e.id];
            }    
            else
                sta[e.id][(tau + j)/period][day] += 1 / cap[e.id];
            stat_Q[tau + j][day] += 1;
            maxcong = max(maxcong, sta[e.id][(tau + j)/period][day]);
        }
        v = e.to;
    }
    return;
}
//compute the standard deviation
double stdab(cand a, cand b){
    double cov = 0;
    for (int i = 0; i < MAX_DAY;i++){
        cov += (a.rec[i] - a.ha) * (b.rec[i] - b.ha);
    }
    cov = cov / (MAX_DAY - 1);
    return sqrt(a.std * a.std + b.std * b.std - 2 * cov);
}

FILE *fp_network, *fp_candidate;
//process different datasets
void preprocess(string query){
    fscanf(fp_network,"%d%d", &V, &M);
    int a, b;
    double c;
    for (int i = 0; i < M; i++) {
        fscanf(fp_network,"%d%d%lf", &a, &b, &c);
        G1[a].push_back(edge(b, i, c/60));
        G2[b].push_back(edge(a, i, c/60));
        cap[i] = c*40000/3600/50*2;
        if(cap[i]<2)
            cap[i] = 100;
    }
    for (int i = 0; i < MAX_DAY;i++){
        printf("processing %d \n", i);
        string ofilename;
        if(query.find("Didi")!= string::npos)
            ofilename = query+to_string(i+20161001) + string(" 0 1440_nid");
        else if(query.find("Beijing")!= string::npos)
            ofilename = query+to_string(i+20170401) + string(" 0 1440_nid");
        else
            ofilename = query+to_string(i+20130801) + string(" 0 1440_nid");
        FILE *fp = fopen(ofilename.c_str(), "r");
        int nq;
        Lambda = 1;
        fscanf(fp, "%d%lf", &nq, &detour);
        memset(load, 0, sizeof(load));
        detour = det;
        for (int j = 0; j < nq;j++){
            int t, s, d;
            fscanf(fp, "%d%d%d", &t, &s, &d);
            int back = 60;
            if(query.find("Didi")!= string::npos)
                back=60;
            else if(query.find("Beijing")!= string::npos)
                back=30;
            else
                back=60;
            if(t<t_lb-back)////////////////////////////////////
                continue;
            if(t>=t_ub)
                break;
            sor_stat(i, t, s, d);
        }
    }
    map<int, map<int,double[MAX_DAY]> >::iterator muliter;
    map<int,double[MAX_DAY]>::iterator iter;
    for (muliter = sta.begin(); muliter != sta.end();muliter++)
    {
        int eid = muliter->first;
        for (iter = muliter->second.begin(); iter != muliter->second.end();iter++){
            int tau = iter->first;
            sta_alg2.push_back(cand(eid, tau, iter->second));
        }
    }
    sort(sta_alg2.begin(), sta_alg2.end(), Cmp);
    double lb = 0;
    set<PI> alg3, alg4;
    for (int i = 0; i < sta_alg2.size();i++){
        lb = max(lb, sta_alg2[i].ha - sta_alg2[i].radius);
        if(sta_alg2[i].ha + sta_alg2[i].radius < lb)
            break;
        e_alg3[sta_alg2[i].t][sta_alg2[i].eid] = i;
        t_alg3[sta_alg2[i].eid][sta_alg2[i].t] = i;
        alg3.insert(PI(sta_alg2[i].eid, sta_alg2[i].t));
    }
    for (int l = 0; l < MAX_T; l++){
        if(e_alg3[l].size()==0)
            continue;
        for (int i = 0; i < V; i++)
        {
            for (int j = 0; j < G1[i].size(); j++)
            {
                if (e_alg3[l].find(G1[i][j].id) == e_alg3[l].end())
                    continue;
                for (int k = j + 1; k < G1[i].size(); k++)
                {
                    if (e_alg3[l].find(G1[i][k].id) == e_alg3[l].end())
                        continue;
                    int ai = e_alg3[l][G1[i][j].id], bi = e_alg3[l][G1[i][k].id];
                    cand a = sta_alg2[ai], b = sta_alg2[bi];
                    if ((a.ha - b.ha + stdab(a, b)) < theta)
                    {
                        if (alg3.count(PI(a.eid,a.t))==1){
                            alg3.insert(PI(b.eid, b.t));
                            et2ha[b.eid][b.t] = et2ha[a.eid][a.t];
                        }
                        if (alg3.count(PI(b.eid,b.t))==1){
                            alg3.insert(PI(a.eid, a.t));
                            et2ha[a.eid][a.t] = et2ha[b.eid][b.t];
                        }
                    }
                }
            }
        }
    }
    map<int, map<int, int>>::iterator it;
    for (it = t_alg3.begin(); it != t_alg3.end();it++){
        for (int tau1 = t_lb; tau1 < t_ub; tau1++){
            int tau2 = tau1 + 1;
            if (it->second.find(tau1) == it->second.end()||it->second.find(tau2) == it->second.end())
                continue;
            int eid = it->first;
            int ai = it->second[tau1], bi = it->second[tau2];
            cand a = sta_alg2[ai], b = sta_alg2[bi];
            if ((a.ha-b.ha+stdab(a,b))<theta){
                if (alg3.count(PI(eid,a.t))==1){
                    alg3.insert(PI(eid, b.t));
                    et2ha[eid][b.t] = et2ha[eid][a.t];
                }
                if (alg3.count(PI(eid,b.t))==1){
                    alg3.insert(PI(eid, a.t));
                    et2ha[eid][a.t] = et2ha[eid][b.t];
                }
            }
        }
    }
    set<PI>::iterator itset;
    for (itset = alg3.begin(); itset != alg3.end();itset++){
        int tau = itset->second, eid = itset->first;
        if (tau < t_lb/period || tau >= t_ub/period)
            continue;
        alg4.insert(*itset);
    }
    fprintf(fp_candidate, "%d %d %d\n", alg4.size() * period, t_lb, t_ub - t_lb);
    for (int i = t_lb; i < t_ub;i++){
        double sum = 0, R_Q = 0;
        for (int j = 0; j < MAX_DAY;j++){
            sum += stat_Q[i][j];
            R_Q = max(R_Q, (double)stat_Q[i][j]);
        }    
        double radius_Q = sqrt(R_Q * R_Q * log(2 / delta) / (2 * MAX_DAY));
        fprintf(fp_candidate, "%f %f\n", sum/MAX_DAY, radius_Q);
    }
    for (itset = alg4.begin(); itset != alg4.end();itset++){
        int tau = itset->second, eid = itset->first;
        if (tau < t_lb/period || tau >= t_ub/period)
            continue;
        for (int j = tau * period; j < (tau + 1) * period; j++)
            fprintf(fp_candidate, "%d %d %f\n", eid, j, et2ha[eid][tau]);
    }
}

int main(int argc, char * argv[]){
    string query;
    if (argc > 1)
        fp_network = fopen(argv[1], "r");
    else
        fp_network = fopen("../data/network_NYC.txt", "r");
    if (argc > 2)
        query = string(argv[2]);
    else
        query = string("../history/NYC/");
    if (argc > 3)
        t_lb = stod(string(argv[3]));
    else
        t_lb = 1020;
    if (argc > 4)
        t_ub = stod(string(argv[4]));
    else
        t_ub = t_lb + 60;
    if (argc > 5)
        fp_candidate = fopen(argv[5], "w");
    else
        fp_candidate = fopen("candidate.txt","w");
    if (argc > 6)
        det = stod(string(argv[6]));
    else
        det = 0.05;
    clock_t proc_s = clock();
    preprocess(query);
    clock_t proc_e = clock();
    double proc_d = proc_e - proc_s;
    #ifdef WATCH_MEM
        watchSolutionOnce(getpid(), usedMemory);
    #endif
    #ifdef WATCH_MEM
        printf("%f %f\n", (double)proc_d / CLOCKS_PER_SEC, usedMemory / 1024.0);
    #else
        printf("%f\n", (double)proc_d / CLOCKS_PER_SEC);
    #endif

    return 0;
}
