#include "Graph.h"


int main() {
    Graph<int> graph(13);


    graph.addEdge(0, 5);
    graph.addEdge(4, 3);
    graph.addEdge(0, 1);
    graph.addEdge(9, 12);
    graph.addEdge(6, 4);
    graph.addEdge(5, 4);
    graph.addEdge(0, 2);
    graph.addEdge(11, 12);
    graph.addEdge(9, 10);
    graph.addEdge(0, 6);
    graph.addEdge(7, 8);
    graph.addEdge(9, 11);
    graph.addEdge(5, 3);
    graph.addEdge(6, 10);
    graph.addEdge(3, 7);

//    graph.print();
//    graph.printAllGraphData();
//    cout << endl;
//    cout << "BFS: " << endl;
//    graph.bfs(2);
//    cout << endl;
//    cout << "DFS:" << endl;
//    graph.dfs(3);
//    cout << endl;
//    vector<vector<int>> shortestPaths = graph.shortestPaths(0);
//    cout << "Shortest paths from 0: " << endl;
//    for (int i = 0; i < shortestPaths.size(); i++) {
//        cout << "From 0 to " << i << ": ";
//        for (int j = 0; j < shortestPaths[i].size(); j++) {
//            cout << shortestPaths[i][j] << " ";
//        }
//        cout << endl;
//    }
//    vector<int> shortestPath = graph.shortestPath(0, 12,true);
//    cout << "Shortest path from 0 to 12: ";
//    for (int i = 0; i < shortestPath.size(); i++) {
//        cout << shortestPath[i] << " ";
//    }
//    // does graph have cycle?
//    cout << "Does graph have cycle? \n" << graph.hasCycle() << endl;
//
//    // check if graph is connected
//    cout << "Is graph connected? \n" << graph.isConnected() << endl;
//    // check how many components are in graph
//    cout << "How many components are in graph? \n" << graph.getNumberOfComponents() << endl;
//    // save all the nodes of each component in a vector of vectors
//    vector<vector<int>> components = graph.getComponents();
//    for (int i = 0; i < components.size(); i++) {
//        cout << "Component " << i << ": ";
//        for (int j = 0; j < components[i].size(); j++) {
//            cout << components[i][j] << " ";
//        }
//        cout << endl;
//    }

    graph.iteratorTest();



    return 0;
}