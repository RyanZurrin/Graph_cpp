#include "Graph.h"


int main() {
//    pair<int, string> p1(1, "A");
//    pair<int, string> p2(2, "B");
//    pair<int, string> p3(3, "C");
//    pair<int, string> p4(4, "D");
//    pair<int, string> p5(5, "E");
//    pair<int, string> p6(6, "F");
//    pair<int, string> p7(7, "G");
//    pair<int, string> p8(8, "H");
//    pair<int, string> p9(9, "I");
//    pair<int, string> p10(10, "J");
//
//
//    Graph<string> graph;
//    graph.addVertex(p1);
//    graph.addVertex(p2);
//    graph.addVertex(p3);
//    graph.addVertex(p4);
//    graph.addVertex(p5);
//    graph.addVertex(p6);
//    graph.addVertex(p7);
//    graph.addVertex(p8);
//    graph.addVertex(p9);
//    graph.addVertex(p10);
//
//
//    graph.addEdge(p1, p2, 2.5);
//    graph.addEdge(p7, p3, 1.0);
//    graph.addEdge(p2, p3, 1.0);
//    graph.addEdge(p3, p4, 1.0);
//    graph.addEdge(p4, p5, 1.0);
//    graph.addEdge(p5, p6, 1.0);
//    graph.addEdge(p6, p7, 1.0);
//    graph.addEdge(p7, p8, 1.0);
//    graph.addEdge(p8, p9, 1.0);
//    graph.addEdge(p9, p10, 1.0);
//    graph.addEdge(p10, p1, 1.0);
//    graph.addEdge(p4, p7, 1.0);
//    graph.addEdge(p8, p1, 1.0);
//    graph.addEdge(p8, p9, 1.0);
//    graph.addEdge(p6, p2, 1.0);
//    graph.addEdge(p10, p1, 1.0);
//
//
//    graph.printNodeData();

    Graph<int> graph(13, true);
//    graph.addEdge(0, 2, 1);
//    graph.addEdge(0, 3, 5);
//    graph.addEdge(0, 4, 2);
//    graph.addEdge(1, 2, 3);
//    graph.addEdge(1, 3, 4);
//    graph.addEdge(1, 4, 6);
//
//
    graph.addEdge(0, 5, 4);
    graph.addEdge(4, 3, 8);
    graph.addEdge(0, 1, 8);
    graph.addEdge(9, 7, 7);
    graph.addEdge(6, 4, 1);
    graph.addEdge(5, 4, 5);
    graph.addEdge(0, 2, 1);
    graph.addEdge(11, 12, 1);
    graph.addEdge(9, 10, 2);
    graph.addEdge(0, 6, 2);
    graph.addEdge(7, 8, 7);
    graph.addEdge(9, 11, 6);
    graph.addEdge(5, 3, 2);
    graph.addEdge(6, 10, 4);
    graph.addEdge(3, 7, 9);
    graph.addEdge(8, 12, 5);
    graph.addEdge(0, 3, 8);
    graph.addEdge(3, 2, 11);
    graph.addEdge(3, 1, 4);
    graph.addEdge(1, 2, 7);
    graph.addEdge(2, 1, 6);
    graph.addEdge(2, 5, 3);
    graph.addEdge(2, 6, 2);
    graph.addEdge(3, 6, 21);
    graph.addEdge(5, 7, 4);
    graph.addEdge(6, 5, 19);
    graph.addEdge(1, 6, 8);
    graph.addEdge(12, 9, 14);
    graph.addEdge(10, 11, 10);
    graph.addEdge(11, 1, 11);
    graph.addEdge(4,11, 6);
    graph.addEdge(4,12, 7);
    graph.addEdge(4,2, 8);

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
    vector<int> shortestPath = graph.shortestPath(0, 12,true);
    cout << "Shortest path from 0 to 12: ";
    for (int i = 0; i < shortestPath.size(); i++) {
        cout << shortestPath[i] << " ";
    }
    // does graph have cycle?
    cout << "Does graph have cycle? \n" << graph.cycle() << endl;

    // check if graph is connected
    cout << "Is graph connected? \n" << graph.connected() << endl;
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
    cout << "Is graph bipartite? \n" << graph.bipartite() << endl;

    //  graph.iteratorTest();
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
    auto testnode =  graph.getVertex(3);
    cout << testnode->getData() << endl;
    // test djikstra
    cout << "Testing djikstra" << endl;
    vector<pair<int, int>> shortestPaths2 = graph.dijkstra(0, 12, true);
    for (int i = 0; i < shortestPaths2.size(); i++) {
        cout << "From 0 to " << i << ": " << shortestPaths2[i].first << " " <<
             shortestPaths2[i].second << endl;
    }
    double shortestPath2 = graph.dijkstra(1, 12);
    cout << "Shortest path from 1to 12: " << shortestPath2 << endl;
    //checking pathLength between 0 and 7
    cout << "Testing pathLength between 0 and 7" << endl;
    cout << graph.pathLength(0, 7) << endl;
    // getting degree of all nodes
    cout << "Testing degree" << endl;
    cout << graph.getDegree(0) << endl;
    cout << graph.getDegree(1) << endl;
    cout << graph.getDegree(2) << endl;
    cout << graph.getDegree(3) << endl;
    cout << graph.getDegree(4) << endl;
    cout << graph.getDegree(5) << endl;
    cout << graph.getDegree(6) << endl;
    cout << graph.getDegree(7) << endl;
    cout << graph.getDegree(8) << endl;
    cout << graph.getDegree(9) << endl;
    cout << graph.getDegree(10) << endl;
    cout << graph.getDegree(11) << endl;
    cout << graph.getDegree(12) << endl;
    // testing getNeighbors
    cout << "Testing getNeighbors" << endl;
    vector<int> neighbors;
    for (int i = 0; i < graph.getV(); ++i) {
        neighbors = graph.getNeighbors(i);
        cout << "Neighbors of " << i << ": ";
        for (int j = 0; j < neighbors.size(); ++j) {
            cout << neighbors[j] << " ";
        }
        cout << endl;
    }
    // set the data in each node to 2+its index
    cout << "Testing setData" << endl;
    for (int i = 0; i < graph.getV(); ++i) {
        graph.setData(2 + i, i);
    }

    // set all the weights to 5
    cout << "Testing setWeight" << endl;
    for (int i = 0; i <= graph.getV(); ++i) {
        for (int j = 0; j <= graph.getV(); ++j) {
            graph.setWeight(i, j, (5.0+i)/(j+1.0));
        }
    }

    graph.printNodeData();
    graph.printAllGraphData();
    // dijkstra
    cout << "Testing dijkstra" << endl;
    double sp2 = graph.dijkstra(1, 12);
    cout << "Shortest path from 1to 12: " << sp2 << endl;
    vector<pair<int, int>> sp = graph.dijkstra(0, 3, true);
    for (int i = 0; i < sp.size(); i++) {
        cout << "From 0 to " << i << ": " << sp[i].first << " " <<
             sp[i].second << endl;
    }



    return 0;
}