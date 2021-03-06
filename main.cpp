#include "Graph.h"
#include <iostream>
#include <chrono>
#include "RunTImer.h"
using namespace std;
#include  "DynamicProblems.h"
int main() {
    RunTimer timer(SECONDS);
    timer.start();
    Graph<int> graph(false);
    graph.addVertex(0);
    graph.addVertex(1);
    graph.addVertex(2);
    graph.addVertex(3);
    graph.addVertex(4);
    graph.addVertex(5);
    graph.addVertex(6);
    graph.addVertex(7);

    graph.addEdge(0, 2);
    graph.addEdge(0, 1);
    graph.addEdge(1, 6);
    graph.addEdge(1, 3);
    graph.addEdge(1, 4);
    graph.addEdge(4, 5);


//    graph.addEdge(0, 1, 3);
//    graph.addEdge(1, 4, 4);
//    graph.addEdge(4, 0, 2);
//    graph.addEdge(1, 5, 1);
//    graph.addEdge(1, 2, 6);
//    graph.addEdge(4, 5, 5);
//    graph.addEdge(2, 6, 2);
//    graph.addEdge(2, 3, 4);
//    graph.addEdge(3, 2, 3);
//    graph.addEdge(3, 7, 6);
//    graph.addEdge(7, 3, 3);
//    graph.addEdge(7, 6, 8);
//    graph.addEdge(5, 6, 4);
//    graph.addEdge(6, 5, 5);


   // graph.addVertex(5);
//    graph.addEdge(0, 1);
//    graph.addEdge(0, 2);
//    graph.addEdge(1, 2);
//    graph.addEdge(2, 3);
//    graph.addEdge(3, 0);
//    graph.addEdge(5, 7, 6);
//    graph.addEdge(7, 5, 4);
//    graph.addEdge(0, 2, 5);
//    graph.addEdge(0, 3, 3);
//    graph.addEdge(3, 4, 2);
//    graph.addEdge(4, 5, 6);
//    graph.addEdge(5, 6, 4);
//    graph.addEdge(6, 7, 3);
//    graph.addEdge(7, 4, 8);
//    graph.addEdge(4, 3, 1);
//    graph.addEdge(3, 2, 3);
//    graph.addEdge(2, 1, 8);
//    graph.addEdge(1, 0, 3);
//    graph.addEdge(0, 1, 11);
//    graph.addEdge(1, 2, 8);
//    graph.addEdge(2, 3, 13);
//    graph.addEdge(4, 5);
//    //graph.addEdge(3, 4);
//    graph.addEdge(5, 4);

    graph.print();
    vector<int> path = graph.topologicalSort();
    for (int i = 0; i < path.size(); i++) {
        cout << path[i] << " ";
    }
    cout << endl;
    // find the strongly connected components
    vector<vector<int>> scc = graph.stronglyConnectedComponents();
    for (int i = 0; i < scc.size(); i++) {
        cout << "SCC " << i << ": ";
        for (int j = 0; j < scc[i].size(); j++) {
            cout << scc[i][j] << " ";
        }
        cout << endl;
    }
    // get the adjacency list representation of the graph
    vector<vector<int>> adjList = graph.getAdjMatrix();
    for (int i = 0; i < adjList.size(); i++) {
        cout << "Vertex " << i << ": ";
        for (int j = 0; j < adjList[i].size(); j++) {
            cout << adjList[i][j] << " ";
        }
        cout << endl;
    }
    // dest djikstra from o to 7
    vector<pair<int, int>> path2 = graph.dijkstra(0, 7, true);
    // is the graph a tree?
    cout << "Is the graph a tree? " << graph.isTree() << endl;
    // use bellman-ford to find the shortest path from vertex 0 to vertex 7
    vector<pair<int, int>>  shortestPath = graph.bellmanFord(0, true);
    // print the bridges
    vector<pair<int, int>> bridges = graph.getBridges();
    for (int i = 0; i < bridges.size(); i++) {
        cout << "Bridge: " << bridges[i].first << " " << bridges[i].second << endl;
    }
    // get articulation points
    vector<int> articulationPoints = graph.getArticulationPoints();
    for (int i = 0; i < articulationPoints.size(); i++) {
        cout << "Articulation point: " << articulationPoints[i] << endl;
    }

    timer.stop();
    timer.display();




    return EXIT_SUCCESS;
}