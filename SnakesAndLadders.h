//
// Created by Ryan.Zurrin001 on 1/9/2022.
//

#ifndef GRAPH_CPP_SNAKESANDLADDERS_H
#define GRAPH_CPP_SNAKESANDLADDERS_H
#include <bits/stdc++.h>
using namespace std;


template <typename T>
class SnL_Graph {
private:
    int V{};
    int E{};
    list<T> *adj;
public:
    explicit SnL_Graph(int V);
    void addEdge(int v, int w, bool undirected = false);
    int minimumCost(int s, int t);
    void printAdjList();
};

template<typename T>
SnL_Graph<T>::SnL_Graph(int V) {
    this->V = V;
    this->E = 0;
    adj = new list<T>[V];
}

template<typename T>
void SnL_Graph<T>::addEdge(int v, int w, bool undirected) {
    adj[v].push_back(w);
    if (undirected) {
        adj[w].push_back(v);
    }
    E++;
}

template<typename T>
int SnL_Graph<T>::minimumCost(int s, int t) {
    queue<int> q;
    vector<bool> visited(V, false);
    vector<int> dist(V, 0);
    q.push(s);
    visited[s] = true;
    dist[s] = 0;
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (auto v : adj[u]) {
            if (!visited[v]) {
                q.push(v);
                visited[v] = true;
                dist[v] = dist[u] + 1;
            }
        }
    }
    return dist[t];
}

template<typename T>
void SnL_Graph<T>::printAdjList() {
    for (int i = 0; i < V; i++) {
        cout << i << "-->";
        for (auto it : adj[i]) {
            cout << it << ", ";
        }
        cout << endl;
    }
}

int min_dice_throws(int n,vector<pair<int,int> > snakes, vector<pair<int,int> > ladders){
    vector<int> dp(n+1,0);
    for(auto s : snakes) {
        int begin = s.first;
        int end = s.second;
        dp[begin] = end - begin;
    }
    for(auto l : ladders) {
        int begin = l.first;
        int end = l.second;
        dp[begin] = end - begin;
    }
    SnL_Graph<int> g(n+1);
    for(int i = 1; i < n; i++) {
        for(int j = 1; j <= 6; j++) {
            int next = i + j;
            next += dp[next];
            if(next <= n) {
                g.addEdge(i, next);
            }
        }
    }
    g.printAdjList();
    return g.minimumCost(1, n);
}
#endif //GRAPH_CPP_SNAKESANDLADDERS_H
