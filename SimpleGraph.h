//
// Created by Ryan.Zurrin001 on 1/9/2022.
//

#ifndef GRAPH_CPP_SIMPLEGRAPH_H
#define GRAPH_CPP_SIMPLEGRAPH_H
#include <bits/stdc++.h>
using namespace std;


template <typename T>
class SimpleGraph {
private:
    int V{};
    int E{};
    list<T> *adj;
public:
    explicit SimpleGraph(int V);
    void addEdge(int v, int w, bool undirected = true);
    void printAdjList();
};

template<typename T>
SimpleGraph<T>::SimpleGraph(int V) {
    this->V = V;
    this->E = 0;
    adj = new list<T>[V];
}

template<typename T>
void SimpleGraph<T>::addEdge(int v, int w, bool undirected) {
    adj[v].push_back(w);
    if (undirected) {
        adj[w].push_back(v);
    }
    E++;

}

template<typename T>
void SimpleGraph<T>::printAdjList() {
    for (int i = 0; i < V; i++) {
        cout << i << "-->";
        for (auto it : adj[i]) {
            cout << it << ", ";
        }
        cout << endl;
    }
}


#endif //GRAPH_CPP_SIMPLEGRAPH_H
