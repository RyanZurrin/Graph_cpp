#include "Graph.h"
#include "BoardGame.h"


int main() {

    vector<pair<int, int>> edges{
            {0,  1},
            {1,  2},
            {0,  4},
            {3,  6},
            {3,  4},
            {4,  5},
            {5,  3},
            {5,  6}
    };
    bool cycle = contains_cycle(7, edges);
    cout << "Cycle: " << cycle << endl;

    // check if the graph is connected

    Graph<int> graph;
    graph.addVertex(0);
    graph.addVertex(1);
    graph.addVertex(2);
    graph.addVertex(3);
    graph.addVertex(4);
    graph.addVertex(5);
    //graph.addVertex(6);

    graph.addEdge(0, 1, 1);
    graph.addEdge(1, 5, 1);
    graph.addEdge(2, 5, 4);
    //graph.addEdge(3, 4, 1);
    //graph.addEdge(2, 3, 9);
    graph.addEdge(2, 0, 6);
    //graph.addEdge(0, 3, 7);
    //graph.addEdge(3, 2, 2);
    graph.addEdge(3, 4, 3);
    //graph.addEdge(4, 5, 10);
    // graph.addEdge(0, 4, 8);
    //graph.addEdge(3, 6, 5);

    graph.print();
    graph.printAllGraphData();
    cout << endl;
    cout << "BFS: " << endl;
    graph.bfs(0);
    cout << endl;
    cout << "DFS:" << endl;
    graph.dfs(0);
    cout << endl;
    vector<vector<int>> shortestPaths = graph.shortestPaths(0, true);
    vector<int> shortestPath = graph.shortestPath(1, 4,true);
    // does graph have cycle?
    cout << "Does graph have cycle? \n" << graph.hasCycle() << endl;

    // check if graph is connected
    cout << "Is graph connected? \n" << graph.isConnected() << endl;
    vector<pair<int, int>> backedges;
    graph.findBackEdges(backedges);
    cout << "Back edges: " << endl;
    for (auto edge : backedges) {
        cout << edge.first << " " << edge.second << endl;
    }

    // check each node to see if there is a cycle in it
    for (int i = 0; i < graph.getV(); i++) {
        cout << "Does node " << i << " have a cycle? \n" <<
             graph.cycleFromVertex(i) << endl;
    }
    // check the weights between each pair of vertices
    for (int i = 0; i < graph.getV(); i++) {
        for (int j = 0; j < graph.getV(); j++) {
            if(graph.getWeight(i, j) != -1) {
                cout << "Weight between " << i << " and " << j << " is "
                     << graph.getWeight(i, j) << endl;
            }

        }
    }
    // print all the ids of all the vertices
    cout << "All the ids of all the vertices: " << endl;
    for (int i = 0; i < graph.getV(); i++) {
        cout << "id of node "<< i << " is " << graph.getId(i) << endl;
    }

    cout<< "\ndijkstras: \n";
    vector<pair<int,int>> path = graph.dijkstra(0, 4, true);
    cout << "\npath: " << endl;
    for (auto & i : path) {
        cout << i.first << " " << i.second << endl;
    }
    // dijkstra
    cout << "\nDijkstra: " << graph.dijkstra(0, 4) << endl;
    unordered_map<int, GraphNode<int>*> nodes = graph.getNodes();
    for (auto node : nodes) {
        cout << "node id: " << node.first << endl;
        cout << "node data: " << node.second->data << endl;
        cout << "node edges: " << endl;
        for (auto edge : node.second->neighbors) {
            cout << edge.first << " " << edge.second << endl;
        }
    }

    // determine if the graph is bipartite
    cout << "\nIs graph bipartite? " << graph.isBipartite() << endl;


    return 0;
}

/*
 Graph<int> graph;
    graph.addVertex(0);
    graph.addVertex(1);
    graph.addVertex(2);
    graph.addVertex(3);
    graph.addVertex(4);
    graph.addVertex(5);
    graph.addVertex(6);

    graph.addEdge(0, 1, true);
    graph.addEdge(1, 2, true);
    graph.addEdge(2, 3, true);
    graph.addEdge(3, 5, true);
    graph.addEdge(5, 6, true);
    graph.addEdge(4, 5, true);
    graph.addEdge(0, 4, true);
    graph.addEdge(3, 4, true);

    graph.print();
    graph.bfs(1);
    cout << endl;

    vector<vector<int>> shortestPaths = graph.shortestPaths(1, true);
    vector<int> shortestPath = graph.shortestPath(1, 6, true);
    cout << "Shortest Paths: " << endl;

    // shortest path from 1 to 4

    vector<pair<int, int>> snakes{
            {34, 12},
            {32, 30},
            {24, 16},
            {20, 6},
            {17, 4}
    };
    vector<pair<int, int>> ladders{
            {2, 15},
            {5, 7},
            {9, 27},
            {18, 29},
            {25, 35}
    };
    int N = 36;
    int minDiceThrows = min_dice_throws(N, snakes, ladders);
    cout << "\nMinimum number of dice throws required to reach the end of the board is " << minDiceThrows << endl;



 */