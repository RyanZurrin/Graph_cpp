#include "Graph.h"


int main() {
    Graph<int> graph(6, true);
    graph.addEdge(0, 1, 1);
    graph.addEdge(0, 3, 5);
    graph.addEdge(2, 1, 2);
    graph.addEdge(5, 1, 3);
    graph.addEdge(4, 2, 4);
    graph.addEdge(4, 3, 6);
    graph.addEdge(5, 1, 7);
    graph.addEdge(5, 3, 8);


//    graph.addEdge(0, 5, 4);
//    //graph.addEdge(4, 3, 8);
//    graph.addEdge(0, 1, 8);
//    graph.addEdge(9, 7, 7);
//    graph.addEdge(6, 4, 1);
//    //graph.addEdge(5, 4, 5);
//    graph.addEdge(0, 2, 1);
//    //graph.addEdge(11, 12, 1);
//    //graph.addEdge(9, 10, 2);
//    //graph.addEdge(0, 6, 2);
//    graph.addEdge(7, 8, 7);
//    graph.addEdge(9, 11, 6);
//    //graph.addEdge(5, 3, 2);
//    graph.addEdge(6, 10, 4);
//    graph.addEdge(3, 7, 9);
//    graph.addEdge(8, 12, 5);
//    graph.addEdge(0, 3, 8);
//    //graph.addEdge(3, 2, 11);
//    //graph.addEdge(3, 1, 4);
//    //graph.addEdge(1, 2, 7);
//    //graph.addEdge(2, 1, 6);
//    //graph.addEdge(2, 5, 3);
//    graph.addEdge(2, 6, 2);
//    graph.addEdge(3, 6, 21);
//    graph.addEdge(5, 7, 4);
//    graph.addEdge(6, 5, 19);
//    graph.addEdge(1, 6, 8);
//    graph.addEdge(12, 9, 14);
//    graph.addEdge(10, 11, 10);
//    graph.addEdge(11, 1, 11);
//    graph.addEdge(4,11, 6);

    graph.print();
    graph.printAllGraphData();
    cout << endl;
    cout << "BFS: " << endl;
    graph.bfs(0);
    cout << endl;
    cout << "DFS:" << endl;
    graph.dfs(0);
    cout << endl;
    vector<vector<int>> shortestPaths = graph.shortestPaths(0);
    cout << "Shortest paths from 0: " << endl;
    for (int i = 0; i < shortestPaths.size(); i++) {
        cout << "From 0 to " << i << ": ";
        for (int j = 0; j < shortestPaths[i].size(); j++) {
            cout << shortestPaths[i][j] << " ";
        }
        cout << endl;
    }
    vector<int> shortestPath = graph.shortestPath(0, 2,true);
    cout << "Shortest path from 0 to 12: ";
    for (int i = 0; i < shortestPath.size(); i++) {
        cout << shortestPath[i] << " ";
    }
    // does graph have cycle?
    cout << "Does graph have cycle? \n" << graph.hasCycle() << endl;

    // check if graph is connected
    cout << "Is graph connected? \n" << graph.isConnected() << endl;
    // check how many components are in graph
    cout << "How many components are in graph? \n" << graph.getNumberOfComponents() << endl;
    // save all the nodes of each component in a vector of vectors
    vector<vector<int>> components = graph.getComponents();
    for (int i = 0; i < components.size(); i++) {
        cout << "Component " << i << ": ";
        for (int j = 0; j < components[i].size(); j++) {
            cout << components[i][j] << " ";
        }
        cout << endl;
    }
    // check if the graph is bipartite
    cout << "Is graph bipartite? \n" << graph.isBipartite() << endl;

    graph.iteratorTest();
    // test each of the graph's methods
    cout << "Testing each of the graph's methods" << endl;
    cout << "Testing addEdge" << endl;
    graph.addEdge(0, 2);
    graph.print();
    cout << "Testing removeEdge" << endl;
    graph.removeEdge(3, 0);
    graph.print();
    cout << "Testing addVertex" << endl;
    graph.addVertex(4);
    graph.print();
    cout << "Testing removeVertex" << endl;
    graph.removeVertex(4);
    graph.print();
    cout << "Testing getNumberOfVertices" << endl;
    cout << graph.getV() << endl;
    cout << "Testing getNumberOfEdges" << endl;
    cout << graph.getE() << endl;
    cout << "Testing getVertex" << endl;
    GraphNode<int>* testnode =  graph.getVertex(3);
    cout << testnode->getData() << endl;
    // test djikstra
    cout << "Testing djikstra" << endl;
    vector<pair<int, int>> shortestPaths2 = graph.dijkstra(0, 3, true);
    for (int i = 0; i < shortestPaths2.size(); i++) {
        cout << "From 0 to " << i << ": " << shortestPaths2[i].first << " " << shortestPaths2[i].second << endl;
    }
    int shortestPath2 = graph.dijkstra(1, 3);
    cout << "Shortest path from 0 to 12: " << shortestPath2 << endl;
    // checking pathLength between 0 and 7
//    cout << "Testing pathLength between 0 and 7" << endl;
//    cout << graph.pathLength(0, 7) << endl;
//    // getting degree of all nodes
//    cout << "Testing degree" << endl;
//    cout << graph.degree(0) << endl;
//    cout << graph.degree(1) << endl;
//    cout << graph.degree(2) << endl;
//    cout << graph.degree(3) << endl;
//    cout << graph.degree(4) << endl;
//    cout << graph.degree(5) << endl;
//    cout << graph.degree(6) << endl;
//    cout << graph.degree(7) << endl;
//    cout << graph.degree(8) << endl;
//    cout << graph.degree(9) << endl;
//    cout << graph.degree(10) << endl;
//    cout << graph.degree(11) << endl;
//    cout << graph.degree(12) << endl;
    // testing getNeighbors
    cout << "Testing getNeighbors" << endl;
    vector<int> neighbors = graph.getNeighbors(0);
    for (int i = 0; i < neighbors.size(); i++) {
        cout << neighbors[i] << " ";
    }
    graph.printAllGraphData();



    return 0;
}