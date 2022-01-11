//
// Created by Ryan.Zurrin001 on 1/9/2022.
//

#ifndef GRAPH_CPP_WEIGHTEDGRAPH_H
#define GRAPH_CPP_WEIGHTEDGRAPH_H
#include<iostream>
#include<list>
#include<set>
using namespace std;

class Graph{

    int V;
    list<pair<int,int> > *l;

public:
    explicit Graph(int v){
        V = v;
        l = new list<pair<int,int> >[V];
    }

    void addEdge(int u,int v,int wt,bool undirected = true){

        l[u].emplace_back(wt,v);
        if(undirected){
            l[v].emplace_back(wt,u);
        }
    }

};
#endif //GRAPH_CPP_WEIGHTEDGRAPH_H
