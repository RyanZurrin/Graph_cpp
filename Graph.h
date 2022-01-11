#pragma ide diagnostic ignored "bugprone-branch-clone"
#pragma ide diagnostic ignored "readability-use-anyofallof"
#pragma ide diagnostic ignored "misc-no-recursion"
#ifndef GRAPH_CPP_GRAPH_H
#define GRAPH_CPP_GRAPH_H
#include <bits/stdc++.h>
#include "UnionFind.h"
using namespace std;
// ___________________________begin graph node class____________________________
/**
 * @brief The GraphNode class
 * @tparam T  the type of the node
 */
template <class T>
class GraphNode {
public:
    T data;
    list<pair<T, double>> neighbors;
    T* parent;
    bool visited{};
    double distance{};
    GraphNode() : data(0), parent(nullptr), visited(false), distance(0) {}
    explicit GraphNode(T data) {
        this->data = data;
        this->visited = false;
        this->distance = 0;
        this->parent = nullptr;
    }
    GraphNode(T data, T parent, double distance, bool visited) {
        this->data = data;
        this->parent = parent;
        this->distance = distance;
        this->visited = visited;
    }
    GraphNode(T data, T parent, double distance) {
        this->data = data;
        this->parent = parent;
        this->distance = distance;
        this->visited = false;
    }
    [[maybe_unused]] void setParent(T parent_) { this->parent = parent_; }
    [[maybe_unused]] void setDistance(double distance_) { this->distance = distance_; }
    [[maybe_unused]] void setVisited(bool visited_) { this->visited = visited_; }

}; //___________________________end graph node class____________________________
/**
 * @brief a comparator for the priority queue that compares the weights of two
 *      nodes
 */
template<class T>
class [[maybe_unused]] GraphNodeComparator {
public:
    // compare the weights of two nodes
    bool operator()(const pair<T, double> &a, const pair<T, double> &b) {
        return a.second > b.second;
    }
}; // ___________________end graph node comparator class________________________
// _______________________________begin graph class_____________________________
template <class T>
class [[maybe_unused]] Graph {
    unordered_map<T, GraphNode<T>*> nodes;

    int V{};
    int E{};
    bool isDirected{};
    [[maybe_unused]] void bfsUtil(T src,
                                  [[maybe_unused]] const bool* visited,
                                  [[maybe_unused]] const bool* recur);

    [[maybe_unused]] void dfsUtil([[maybe_unused]] T src,
                                  [[maybe_unused]] bool* visited,
                                  [[maybe_unused]] bool* recur);

    [[maybe_unused]] void dfsUtil(T src,
                                  [[maybe_unused]] bool* visited);

    [[maybe_unused]] bool hasCycleUtil([[maybe_unused]] T node,
                                       bool* visited,
                                       [[maybe_unused]] bool* recur,
                                       [[maybe_unused]] T parent);

    [[maybe_unused]] bool isConnectedUtil([[maybe_unused]] int i,
                                          [[maybe_unused]] bool* visited);
    [[maybe_unused]] void traverse([[maybe_unused]] T u, bool* visited);
    [[maybe_unused]] bool cycleFromVertexUtil(T node, [[maybe_unused]] bool* visited);
    [[maybe_unused]] void findBackEdgesUtil([[maybe_unused]] T node,
                                            [[maybe_unused]] bool* visited,
                                            [[maybe_unused]] bool* recur,
                                            [[maybe_unused]] vector<pair<T,T>>& backEdges);
    [[maybe_unused]] bool isBipartiteUtil(T node, bool* visited, int* color);
    [[maybe_unused]] int encode(int i, int j, int rows, int cols, bool rowMajor = true);

public:
    [[maybe_unused]] explicit Graph(bool isDirected_ = false);
    [[maybe_unused]] explicit Graph(vector<vector<T>> grid,
                                    int numbVertices,
                                    bool isDirected = false);
    [[maybe_unused]] explicit Graph(int vertices, bool isDirected = false);
    [[maybe_unused]] explicit Graph([[maybe_unused]] vector<T> nodes,
                                    [[maybe_unused]] bool isDirected_ = false);
    [[maybe_unused]] explicit Graph(const string& fileName,
                                    char delimiter = ' ',
                                    bool isDirected_ = false);
    [[maybe_unused]] Graph(T* nodes, int size, bool isDirected_ = false);
    [[maybe_unused]] void addVertex(T i);
    [[maybe_unused]] void addEdge(T node1, T node2);
    [[maybe_unused]] void addEdge(T node1, T node2, double weight);
    [[maybe_unused]] void removeEdge(T node1, T node2);
    [[maybe_unused]] void removeNode(T node);
    [[maybe_unused]] void bfs([[maybe_unused]] T src);
    [[maybe_unused]] void dfs(T src);
    [[maybe_unused]] void findBackEdges(vector<pair<T,T>>& backEdges);
    [[maybe_unused]] void print();
    [[maybe_unused]] void printAllGraphData();
    [[maybe_unused]] bool isConnected();
    [[maybe_unused]] bool hasCycle();
    [[maybe_unused]] bool cycleFromVertex(T node);
    [[maybe_unused]] bool isBipartite();
    [[maybe_unused]] bool findEdge(T node1, T node2);
    [[maybe_unused]] [[nodiscard]] bool directed() const;
    [[maybe_unused]] [[nodiscard]] int getV() const;
    [[maybe_unused]] [[nodiscard]] int getE() const;
    [[maybe_unused]] double getWeight(T node1, T node2);
    [[maybe_unused]] T getId(T node);
    [[maybe_unused]] int minDistance(const int *pInt, const bool *pBoolean);
    [[maybe_unused]] [[maybe_unused]] [[nodiscard]] unordered_map<T, GraphNode<T>*> getNodes() const;
    [[maybe_unused]] vector<pair<T,T>> dijkstra(T src, T dest, bool print);
    [[maybe_unused]] [[maybe_unused]] T dijkstra(T src, T dest);
    [[maybe_unused]] vector<T> shortestPath(T src, T dest, bool print = false);
    [[maybe_unused]] vector<vector<T>> shortestPaths(T src, bool print = false);

}; // ____________________________end graph class_______________________________
// ____________________begin graph class function definitions___________________
// *****************************************************************************
/**
 * @brief Graph constructor
 * @param isDirected_ : default is false, set to true if graph is directed
 */
template <class T>
Graph<T>::Graph(bool isDirected_) {
    isDirected = isDirected_;
} // ____________________________end graph constructor__________________________

template<class T>
Graph<T>::Graph(vector<vector<T>> grid, int numbVertices, bool isDirected) {
    // create graph with grid
    int rows = grid.size();
    int cols = grid[0].size();
    this->isDirected = isDirected;
    // add numVertices vertices to graph
    for (int i = 0; i <= numbVertices-1; i++) {
        this->addVertex(i);
    }
    // build  arrays dx and dy to store the directions of the 4-connected neighbors
    int dx[4] = {-1, 0, 1, 0};
    int dy[4] = {0, 1, 0, -1};
    // each pos in grid has 4 neighbors to N E S W if in bounds and the weight
    // is the value at that pos. use the encode to add neighbors to graph and weight
    // to the value of the pos at the N E S W
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < 4; k++) {
                int x = i + dx[k];
                int y = j + dy[k];
                int encodedIJ = encode(i, j, rows, cols);
                int encodedXY = encode(x, y, rows, cols);
                if (x >= 0 && x < rows && y >= 0 && y < cols) {
                    // use findEdge to make sure we don't add a duplicate edge
                    if (!findEdge(encodedIJ, encodedXY)) {
                        this->addEdge(encodedIJ, encodedXY, grid[x][y]);
                    }
                }
            }
        }
    }
}

/**
 * @brief Graph constructor with number of vertices which will be set from
 *  zero to number of vertices - 1
 * @param vertices  : number of vertices
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>[[maybe_unused]]
Graph<T>::Graph(int vertices, bool isDirected_) {
    isDirected = isDirected_;
    for (int i = 0; i < vertices; i++) {
        addVertex(i);
    }
} // ____________________________end graph constructor__________________________
/**
 * @brief Graph constructor with vector of nodes
 * @param nodes  : vector of nodes
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>[[maybe_unused]]
Graph<T>::Graph([[maybe_unused]] vector<T> nodes, [[maybe_unused]] bool isDirected_) {
    isDirected = isDirected_;
    for ([[maybe_unused]] auto node : nodes) {
        this->nodes[node] = new GraphNode<T>(node);
        V++;
    }
} // ____________________________end graph constructor__________________________
/**
 * @brief Graph constructor with a file name and delimiter
 * @param fileName  : file name
 * @param delimiter  : delimiter
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>[[maybe_unused]]
Graph<T>::Graph(const string& fileName,
                [[maybe_unused]] char delimiter,
                bool isDirected_) {
    isDirected = isDirected_;
    ifstream file(fileName);
    if (!file.is_open()) {
        cout << "File not found" << endl;
        return;
    }
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        string node1, node2;
        double weight;
        ss >> node1 >> node2 >> weight;
        addEdge(node1, node2, isDirected);
        if (!isDirected) {
            addEdge(node2, node1, isDirected);
        }
    }
} // ____________________________end graph constructor__________________________
/**
 * @brief Graph constructor with an array of nodes
 * @param nodes  : array of nodes
 * @param size  : size of array
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>[[maybe_unused]]
Graph<T>::Graph(T *nodes, int size, bool isDirected_) {
    isDirected = isDirected_;
    for (int i = 0; i < size; i++) {
        this->nodes[nodes[i]] = new GraphNode<T>(nodes[i]);
        V++;
    }
} // ____________________________end graph constructor__________________________

/**
 * @brief adds a vertex to the graph
 * @param i  : node to be added
 */
template<class T>[[maybe_unused]]
void Graph<T>::addVertex(T i) {
    nodes[i] = new GraphNode<T>(i);
    V++;
} // ____________________________end addVertex__________________________________
/**
 * @brief adds an edge to the graph
 * @param node1  : first node, in directed graph this is the tail
 * @param node2  : second node, in directed graph this is the head\n
 * if directed:  node1 ---> node2
 */
template<class T>
void Graph<T>::addEdge(T node1, T node2) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        throw invalid_argument("One of the nodes does not exist");
    }
    nodes[node1]->neighbors.push_back(make_pair(node2, 0));
    if (!isDirected) {
        nodes[node2]->neighbors.push_back(make_pair(node1, 0));
        E++;
    }
    E++;
} // ___________________________end addEdge_____________________________________
/**
 * @brief adds a weighted edge to the graph
 * @param node1  : first node, in directed graph this is the tail
 * @param node2  : second node, in directed graph this is the head\n
 * @param weight  : weight of the edge\n
 * if directed:  node1--(weight)-->node2
 */
template<class T>[[maybe_unused]]
void Graph<T>::addEdge(T node1, T node2, double weight) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        throw invalid_argument("One of the nodes does not exist");
    }
    nodes[node1]->neighbors.push_back(make_pair(node2, weight));
    if (!isDirected) {
        nodes[node2]->neighbors.push_back(make_pair(node1, weight));
        E++;
    }
    E++;
} // _____________________________end addEdge___________________________________
/**
 * @brief removes an edge from the graph
 * @param node1  : first node
 * @param node2  : second node
 */
template<class T>[[maybe_unused]]
void Graph<T>::removeEdge(T node1, T node2) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        throw invalid_argument("One of the nodes does not exist");
    }
    nodes[node1]->neighbors.remove_if([&](pair<T, double> neighbor) {
        return neighbor.first == node2;
    });
    if (!isDirected) {
        nodes[node2]->neighbors.remove_if([&](pair<T, double> neighbor) {
            return neighbor.first == node1;
        });
    }
    E--;
} // ________________________end removeEdge_____________________________________
/**
 * @brief removes a vertex from the graph
 * @param i  : node to be removed
 */
template<class T>[[maybe_unused]]
void Graph<T>::removeNode([[maybe_unused]] T node) {
    if (nodes.find(node) == nodes.end()) {
        throw invalid_argument("Node does not exist");
    }
    for ([[maybe_unused]] auto neighbor : nodes[node]->neighbors) {
        removeEdge(node, neighbor.first);
    }
    nodes.erase(node);
    V--;
} // ________________________end removeNode_____________________________________
/**
 * @brief breadth first search
 * @param src  : source node to start search
 */
template<class T> [[maybe_unused]]
void Graph<T>::bfs([[maybe_unused]] T src) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    queue<T> q;
    [[maybe_unused]] bool* visited = new bool[V]{false};
    q.push(src);
    visited[src] = true;
    while (!q.empty()) {
        T node = q.front();
        q.pop();
        cout << node << " ";
        for ([[maybe_unused]] auto neighbor : nodes[node]->neighbors) {
            if (!visited[neighbor.first]) {
                q.push(neighbor.first);
                visited[neighbor.first] = true;
            }
        }
    }
    cout << endl;
    delete[] visited;
} // ________________________________end bfs____________________________________
/**
 * @brief depth first search
 * @param src  : source node to start search
 */
template<class T> [[maybe_unused]]
void Graph<T>::dfs([[maybe_unused]] T src) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    [[maybe_unused]] bool* visited = new bool[V]{false};
    dfsUtil(src, visited);
    delete[] visited;
} // _________________________________end dfs___________________________________

template<class T>
void Graph<T>::findBackEdges(vector<pair<T, T>> &backEdges) {
    bool *visited = new bool[V];
    bool *recur = new bool[V];
    for (int i = 0; i < V; i++) {
        visited[i] = false;
        recur[i] = false;
    }
    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            findBackEdgesUtil(i, visited, recur, backEdges);
        }
    }
    delete[] visited;
    delete[] recur;
}

/**
 * @brief prints the graph in adjacency list format
 */
template<class T> [[maybe_unused]]
void Graph<T>::print() {
    for ([[maybe_unused]] auto node : nodes) {
        cout << node.first << ": ";
        for ([[maybe_unused]] auto neighbor : node.second->neighbors) {
            cout << neighbor.first << " ";
        }
        cout << endl;
    }
} // ________________________end print__________________________________________
/**
 * @brief prints the graph in adjacency matrix format with weights
 */
template<class T> [[maybe_unused]]
void Graph<T>::printAllGraphData() {
    cout << "V: " << V << endl;
    cout << "E: " << E << endl;
    // print the weights of the edges as well
    for ([[maybe_unused]] auto node : nodes) {
        cout << node.first << ":";
        for ([[maybe_unused]] auto neighbor : node.second->neighbors) {
            cout << "(" << neighbor.second << ")" << neighbor.first << " ";
        }
        cout << endl;
    }
} // ________________________end printAllGraphData______________________________
/**
 * @brief determines if the graph is connected
 * @return  : true if the graph is connected, false otherwise
 */
template<class T>
bool Graph<T>::isConnected() {
    if (V == 1) {
        return true;
    }
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++) {
        visited[i] = false;
    }
    isConnectedUtil(0, visited);
    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            delete[] visited;
            return false;
        }
    }
    delete[] visited;
    return true;
} // _____________________________end isConnected_______________________________
/**
 * @brief determines if the graph has a cycle
 * @warning : this function does not work for disconnected undirected graphs
 * @return  : true if the graph has a cycle, false otherwise
 */
template<class T>
bool Graph<T>::hasCycle() {
    [[maybe_unused]] bool* visited;
    T parent = -1;
    // set all nodes to unvisited
    visited = new bool[V];
    [[maybe_unused]] bool *recStack = new bool[V];
    for ([[maybe_unused]] int i = 0; i < V; i++) {
        recStack[i] = false;
        visited[i] = false;
    }
    if (isDirected) {
        for ([[maybe_unused]] auto node : nodes) {
            if (hasCycleUtil(node.first, visited, recStack, parent)) {
                delete[] recStack;
                return true;
            }
        }
    } else {
        return hasCycleUtil(0, visited, recStack, parent);
    }
    delete[] recStack;
    return false;
} // _____________________________end hasCycle__________________________________
/**
 * @brief determines if there is a cycle from the given node of a directed graph
 * @warning this function is only valid for directed graphs
 * @param node  : node to start search from
 * @return  : true if there is a cycle, false otherwise
 */
template<class T>
[[maybe_unused]] bool Graph<T>::cycleFromVertex(T node) {
    if (nodes.find(node) == nodes.end()) {
        throw invalid_argument("Node does not exist");
    }
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++) {
        visited[i] = false;
    }
    return cycleFromVertexUtil(node, visited);
} // _____________________________end cycleFromVertex___________________________
/**
 * @brief determines if the graph is bipartite
 * @return  : true if the graph is bipartite, false otherwise
 */
template<class T>
bool Graph<T>::isBipartite() {
    int color[V];
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++) {
        visited[i] = false;
        color[i] = -1;
    }
    for ([[maybe_unused]] auto node : nodes) {
        if (color[node.first] == -1) {
            if (!isBipartiteUtil(node.first, visited, color)) {
                return false;
            }
        }
    }
    return true;
} // _____________________________end isBipartite_______________________________
/**
 * @brief determines if an edge exists between two nodes
 * @param node1  : first node
 * @param node2  : second node
 * @return  : true if an edge exists, false otherwise
 */
template<class T>
bool Graph<T>::findEdge(T node1, T node2) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        throw invalid_argument("Node does not exist");
    }
    // iterate over the neighbors of node1
    for ([[maybe_unused]] auto neighbor : nodes[node1]->neighbors) {
        if (neighbor.first == node2) {
            return true;
        }
    }
    return false;
} // _____________________________end findEdge__________________________________
/**
 * @brief returns the isDirected property of the graph
 * @return  : true if the graph is directed, false otherwise
 */
template<class T>
[[maybe_unused]] bool Graph<T>::directed() const {
    return this->isDirected;
} // ________________________end directed_______________________________________
/**
 * @brief number of vertices in the graph
 * @return  : number of vertices in the graph
 */
template<class T>[[maybe_unused]]
int Graph<T>::getV() const {
    return V;
} // ________________________end getV___________________________________________
/**
 * @brief number of edges in the graph
 * @return  : number of edges in the graph
 */
template<class T>[[maybe_unused]]
int Graph<T>::getE() const {
    return E;
} // ________________________end getE___________________________________________
/**
 * @brief getter function for the nodes of the graph
 * @return  : unordered map of all the nodes in the graph
 */
template<class T>
[[maybe_unused]] unordered_map<T, GraphNode<T> *> Graph<T>::getNodes() const {
    return nodes;
}// _________________________end getNodes_______________________________________

/**
 *  @brief algorithm to find the shortest path between two nodes
 * @param src  : source node
 * @param dest  : destination node
 * @return  : the shortest path between the two nodes
 */
template<class T>[[maybe_unused]]
T Graph<T>::dijkstra([[maybe_unused]] T src,
                     [[maybe_unused]] T dest) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    if (nodes.find(dest) == nodes.end()) {
        throw invalid_argument("Destination node does not exist");
    }
    if (src == dest) {
        throw invalid_argument("Source and destination are the same");
    }
    vector<double> distances(V, INT_MAX);
    set<pair<double, T>> s;
    distances[src] = 0;
    s.insert(make_pair(0, src));
    while (!s.empty()) {
        [[maybe_unused]] auto it = s.begin();
        [[maybe_unused]] auto curNode = it->second;
        [[maybe_unused]] auto distTilNow = it->first;
        s.erase(it);
        for ([[maybe_unused]] auto node : nodes[curNode]->neighbors) {
            [[maybe_unused]] auto neighbor = node.first;
            [[maybe_unused]] auto currEdge = node.second;
            if (distTilNow + currEdge < distances[neighbor]) {
                [[maybe_unused]] auto f = s.find({distances[neighbor], neighbor});
                if (f != s.end()) {
                    s.erase(f);
                }
                distances[neighbor] = distTilNow + currEdge;
                s.insert({distances[neighbor], neighbor});
            }
        }
    }
//    // Single source shortest path to all other nodes
//    for (int i = 0; i < V; ++i) {
//        cout << "Distance of " << i << " is " << distances[i] << endl;
//    }
//    cout  << endl;
    // check if the distance is INT_MAX
    if (distances[dest] == INT_MAX) {
        return -1;
    }
    return distances[dest];
} // ________________________end dijkstra_______________________________________
/**
 * @brief algorithm to find the shortest path between two nodes
 * @param src  : source node
 * @param dest  : destination node
 * @param print  : true if the path should be printed, false otherwise
 * @return  : vector of nodes in the shortest path
 */
template<class T>[[maybe_unused]]
vector<pair<T,T>> Graph<T>::dijkstra([[maybe_unused]] T src,
                                     [[maybe_unused]] T dest,
                                     [[maybe_unused]] bool print) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    if (nodes.find(dest) == nodes.end()) {
        throw invalid_argument("Destination node does not exist");
    }
    if (src == dest) {
        throw invalid_argument("Source and destination are the same");
    }
    vector<double> distances(V, INT_MAX);
    set<pair<double, T>> s;
    distances[src] = 0;
    s.insert(make_pair(0, src));
    vector<pair<T,T>> path;
    while (!s.empty()) {
        [[maybe_unused]] auto it = s.begin();
        [[maybe_unused]] auto curNode = it->second;
        [[maybe_unused]] auto distTilNow = it->first;
        s.erase(it);
        for ([[maybe_unused]] auto node : nodes[curNode]->neighbors) {
            [[maybe_unused]] auto neighbor = node.first;
            [[maybe_unused]] auto currEdge = node.second;
            if (distTilNow + currEdge < distances[neighbor]) {
                [[maybe_unused]] auto f = s.find({distances[neighbor], neighbor});
                if (f != s.end()) {
                    s.erase(f);
                }
                distances[neighbor] = distTilNow + currEdge;
                s.insert({distances[neighbor], neighbor});
                path.push_back({curNode,neighbor});
            }
        }
    }
    if (print) {
        for ([[maybe_unused]] auto p : path) {
            cout << p.first << " -> " << p.second << " = " << distances[p.second] << endl;
        }
    }
    return path;
} // ________________________end dijkstra_______________________________________

/**
 * @brief finds the shortest path between two nodes, if it exists
 * @param src  : source node
 * @param dest  : destination node
 * @param print  : true if the path should be printed, false otherwise
 * @return  : vector of nodes in the shortest path, empty vector if no path exists
 */
template<class T> [[maybe_unused]]
vector<T> Graph<T>::shortestPath([[maybe_unused]] T src,
                                 [[maybe_unused]] T dest,
                                 [[maybe_unused]] bool print) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    if (nodes.find(dest) == nodes.end()) {
        throw invalid_argument("Destination node does not exist");
    }
    if (src == dest) {
        vector<T> path;
        path.push_back(src);
        return path;
    }
    vector<T> path;
    [[maybe_unused]] bool *visited = new bool[V];
    [[maybe_unused]] bool *recStack = new bool[V];
    for ([[maybe_unused]] int i = 0; i < V; i++) {
        visited[i] = false;
        recStack[i] = false;
    }
    queue<T> q;
    [[maybe_unused]] int* dist = new int[V]{0};
    [[maybe_unused]] T* prev = new T[V]{-1};
    q.push(src);
    visited[src] = true;
    dist[src] = 0;
    while (!q.empty()) {
        T u = q.front();
        q.pop();
        for ([[maybe_unused]] auto neighbor : nodes[u]->neighbors) {
            if (!visited[neighbor.first]) {
                visited[neighbor.first] = true;
                dist[neighbor.first] = dist[u] + 1;
                prev[neighbor.first] = u;
                q.push(neighbor.first);
            }
        }
    }
    if (visited[dest]) {
        T u = dest;
        while (u != src) {
            path.push_back(u);
            u = prev[u];
        }
        path.push_back(src);
        reverse(path.begin(), path.end());
        if (print) {
            cout << "Shortest path from " << src << " to " << dest << " is: ";
            for ([[maybe_unused]] auto node : path) {
                cout << node << " ";
            }
            cout << endl;
        }
        delete[] visited;
        delete[] recStack;
        delete[] dist;
        delete[] prev;
        return path;
    } else {
        cout << "No path from " << src << " to " << dest << endl;
        delete[] visited;
        delete[] recStack;
        delete[] dist;
        delete[] prev;
        return path;
    }
} // ________________________end shortestPath___________________________________
/**
 * single source shortest path algorithm (SSSP) to find the shortest path from a
 * source node to all other nodes in a unweighted graph
 * @param src : source node
 * @param print : true if the path should be printed, false otherwise
 * @return : vector of nodes in the shortest path, empty vector if no path exists
 */
template<class T>[[maybe_unused]]
vector<vector<T>> Graph<T>::shortestPaths([[maybe_unused]] T src,
                                          [[maybe_unused]] bool print) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    // determine if the graph is weighted
    [[maybe_unused]] bool weighted = false;
    for ([[maybe_unused]] auto node : nodes) {
        for ([[maybe_unused]] auto neighbor : node.second->neighbors) {
            if (neighbor.second != 0) {
                weighted = true;
                continue;
            }
        }
    }
    if (print && weighted) {
        cout << "***************************************************";
        cout << "\nWARNING:  This search does not incorporate weights."
             << "\nUse Dijkstra's algorithm for valid weighted search." << endl;
        cout << "***************************************************\n";
    }
    vector<vector<T>> paths;
    [[maybe_unused]] bool *visited = new bool[V];
    [[maybe_unused]] bool *recStack = new bool[V];
    for ([[maybe_unused]] int i = 0; i < V; i++) {
        visited[i] = false;
        recStack[i] = false;
    }
    queue<T> q;
    [[maybe_unused]] int* dist = new int[V]{0};
    [[maybe_unused]] T* prev = new T[V]{-1};

    q.push(src);
    visited[src] = true;
    dist[src] = 0;
    while (!q.empty()) {
        T u = q.front();
        q.pop();
        for ([[maybe_unused]] auto neighbor : nodes[u]->neighbors) {
            if (!visited[neighbor.first]) {
                visited[neighbor.first] = true;
                dist[neighbor.first] = dist[u] + 1;
                prev[neighbor.first] = u;
                q.push(neighbor.first);
            }
        }
    }
    for ([[maybe_unused]] int i = 0; i < V; i++) {
        if (visited[i]) {
            vector<T> path;
            T u = i;
            while (u != src) {
                path.push_back(u);
                u = prev[u];
            }
            path.push_back(src);
            reverse(path.begin(), path.end());
            paths.push_back(path);
            if (print) {
                cout << "Shortest path from " << src << " to " << i << " is: ";
                for ([[maybe_unused]] auto node : path) {
                    cout << node << " ";
                }
                cout << endl;
            }
        }
    }
    delete[] visited;
    delete[] recStack;
    delete[] dist;
    delete[] prev;
    return paths;
} // ________________________end shortestPaths__________________________________


// ___________________________bfsUtil___________________________________________
/**
 * @brief private helper function to help with bfs traversals
 * @param src  source node
 * @param visited  vector of visited nodes
 * @param recur  vector representing a stack of nodes to keep track of what node
 *  a node was visited from
 */
template<class T> [[maybe_unused]]
void Graph<T>::bfsUtil(
        [[maybe_unused]] T src, [[maybe_unused]] const bool* visited,
        [[maybe_unused]] const bool* recur) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    queue<T> q;
    q.push(src);
    visited[src] = true;
    recur[src] = true;
    while (!q.empty()) {
        T node = q.front();
        q.pop();
        cout << node << " ";
        for ([[maybe_unused]] auto neighbor : nodes[node]->neighbors) {
            if (!visited[neighbor.first]) {
                q.push(neighbor.first);
                visited[neighbor.first] = true;
                if (!isDirected) {
                    recur[neighbor.first] = true;
                }
            }
        }
    }
} // ________________________end bfsUtil________________________________________

/**
 * @brief private helper function to help with dfs traversals
 * @param src  source node
 * @param visited  vector of visited nodes
 * @param recur vector representing a stack of nodes to keep track of what node
 *  a node was visited from
 */
template<class T>
void Graph<T>::dfsUtil([[maybe_unused]] T src,
                       [[maybe_unused]] bool *visited,
                       [[maybe_unused]] bool *recur) {
    visited[src] = true;
    recur[src] = true;
    for ([[maybe_unused]] auto neighbor : nodes[src]->neighbors) {
        if (!visited[neighbor.first]) {
            dfsUtil(neighbor.first, visited, recur);
        }
    }
    recur[src] = false;
} // ________________________end dfsUtil________________________________________
/**
 * @brief private helper function to help with dfs traversals of undirected
 * graphs
 * @param src  source node
 * @param visited  vector of visited nodes
 */
template<class T> [[maybe_unused]]
void Graph<T>::dfsUtil([[maybe_unused]]  T src,
                       [[maybe_unused]] bool* visited) {
    visited[src] = true;
    for ([[maybe_unused]] auto neighbor : nodes[src]->neighbors) {
        if (!visited[neighbor.first]) {
            dfsUtil(neighbor.first, visited);
        }
    }
    cout << src << " ";
} // ________________________end dfsUtil________________________________________

//___________________________isCyclicUtil_______________________________________
/**
 * @brief private helper function for finding cycles in the graph
 * @param node  node to start the search from
 * @param visited  vector of visited nodes
 * @param recur vector representing a stack of nodes
 * @param parent  parent node of the current node
 * @return  true if a cycle is found, false otherwise
 */
template<class T>[[maybe_unused]]
bool Graph<T>::hasCycleUtil([[maybe_unused]] T node,
                            [[maybe_unused]] bool* visited,
                            [[maybe_unused]] bool* recur,
                            [[maybe_unused]] T parent) {
    if (isDirected) {
        if (!visited[node]) {
            visited[node] = true;
            recur[node] = true;
            // recur all the vertices adjacent to this vertex
            for ([[maybe_unused]] auto neighbor : nodes[node]->neighbors) {
                if (!visited[neighbor.first] &&
                    hasCycleUtil(neighbor.first, visited, recur, parent)) {
                    return true;
                } else if (recur[neighbor.first]) {
                    return true;
                }
            }
        }
        recur[node] = false;
        return false;
    } else {
        visited[node] = true;
        recur[node] = true;
        for ([[maybe_unused]] auto neighbor : nodes[node]->neighbors) {
            if (!visited[neighbor.first]) {
                [[maybe_unused]] bool cycle = hasCycleUtil(neighbor.first, visited, recur, node);
                if (cycle) {
                    return true;
                }
            } else if (neighbor.first != parent) {
                return true;
            }
        }
        return false;
    }
} // ________________________end isCyclicUtil___________________________________
/**
 * @brief private helper function to determine if graph is connected
 * @param i node to start the search from
 * @param visited  vector of visited nodes
 * @param recur  vector representing a stack of nodes
 * @return true if the graph is connected, false otherwise
 */
template<class T>
bool Graph<T>::isConnectedUtil([[maybe_unused]] int i,
                               [[maybe_unused]] bool* visited) {
    visited[i] = true;
    if (!isDirected) { // undirected
        for ([[maybe_unused]] auto neighbor : nodes[i]->neighbors) {
            if (!visited[neighbor.first]) {
                [[maybe_unused]] bool cycle =
                        isConnectedUtil(neighbor.first, visited);
                if (cycle) {
                    return true;
                }
            }
        }
        return false;
    } else { // directed
        traverse(i, visited);
        for ([[maybe_unused]] auto neighbor : nodes[i]->neighbors) {
            if (!visited[neighbor.first]) {
                return false;
            }
        }
        return true;
    }
} // ________________________end isConnectedUtil_______________________________
/**
 * @brief basic traversal of the graph
 * @param u node to look at in the visited vector
 * @param visited  vector of visited nodes
 */
template<class T>
void Graph<T>::traverse([[maybe_unused]] T u, [[maybe_unused]] bool *visited) {
    visited[u] = true;
    for ([[maybe_unused]] auto neighbor : nodes[u]->neighbors) {
        if (!visited[neighbor.first]) {
            traverse(neighbor.first, visited);
        }
    }
} // ________________________end traverse_______________________________________

/**
 * private utility function to determine if graph has a cycle from a given node
 * @param node  node to start the search from
 * @param visited  vector of visited nodes
 * @param parent  parent node of the current node
 */
template<class T>
bool Graph<T>::cycleFromVertexUtil([[maybe_unused]] T node,
                                   [[maybe_unused]] bool *visited) {
    visited[node] = true;
    for ([[maybe_unused]] auto neighbor : nodes[node]->neighbors) {
        if (!visited[neighbor.first]) {
            if (cycleFromVertexUtil(neighbor.first, visited)) {
                return true;
            }
        } else if (neighbor.first != node) {
            return true;
        }
    }
    return false;
}
/**
 * @brief private utility function to determine if graph has a cycle
 * @param node  node to start the search from
 * @param visited  vector of visited nodes
 * @param parent  parent node of the current node
 */
template<class T>
void Graph<T>::findBackEdgesUtil([[maybe_unused]] T node,
                                 [[maybe_unused]] bool *visited,
                                 [[maybe_unused]] bool *recur,
                                 [[maybe_unused]] vector<pair<T, T>> &backEdges) {
    visited[node] = true;
    recur[node] = true;
    for ([[maybe_unused]] auto neighbor : nodes[node]->neighbors) {
        if (!visited[neighbor.first]) {
            findBackEdgesUtil(neighbor.first, visited, recur, backEdges);
        } else if (recur[neighbor.first]) {
            backEdges.push_back(make_pair(node, neighbor.first));
        }
    }
    recur[node] = false;
}

template<class T>[[maybe_unused]]
int Graph<T>::minDistance(const int *pInt, const bool *pBoolean) {
    int min = INT_MAX;
    int min_index;
    for (int i = 0; i < V; i++) {
        if (pBoolean[i] == false && pInt[i] <= min) {
            min = pInt[i];
            min_index = i;
        }
    }
    return min_index;
}

template<class T>[[maybe_unused]]
double Graph<T>::getWeight([[maybe_unused]] T node1, [[maybe_unused]] T node2) {
    // nodes are out of bounds
    if (node1 >= V || node2 >= V) {
        return -1;
    }
    for ([[maybe_unused]] auto neighbor : nodes[node1]->neighbors) {
        if (neighbor.first == node2) {
            return neighbor.second;
        }
    }
    return -1;
}

template<class T>
T Graph<T>::getId(T node) {
    auto it = nodes.find(node);
    if (it != nodes.end()) {
        return it->first;
    }
    return -1;
}

template<class T>
bool Graph<T>::isBipartiteUtil(T node, bool *visited, int *color) {
    color[node] = 1;
    visited[node] = true;
    queue<T> q;
    q.push(node);
    while (!q.empty()) {
        auto u = q.front();
        q.pop();
        for ([[maybe_unused]] auto neighbor : nodes[u]->neighbors) {
            if (!visited[neighbor.first]) {
                if (color[neighbor.first] == color[u]) {
                    return false;
                }
                visited[neighbor.first] = true;
                color[neighbor.first] = !color[u];
                q.push(neighbor.first);
            } else if (color[neighbor.first] == color[u]) {
                return false;
            }
        }
    }
    return true;
}

template<class T>
int Graph<T>::encode(int i, int j, int rows, int cols, bool rowMajor) {
    if (rowMajor) {
        return i * cols + j;
    } else {
        return i + j * rows;
    }
}





// ________________________end cycleFromVertexUtil____________________________

#endif //GRAPH_CPP_GRAPH_H
// ____________________________Static functions_________________________________
/**
 * @brief static function to determine if a graph has a cycle
 * @param V  number of vertices in the graph
 * @param edges  vector of edges in the graph
 * @param directed  default is false, set to true if graph is directed
 * @return true if the graph has a cycle, false otherwise
 */
static bool contains_cycle(int V,
                           const vector<pair<int,int> >& edges,
                           bool directed = false) {
    //Complete this method
    Graph<int> g(directed);
    for (int i = 0; i < V; i++) {
        g.addVertex(i);
    }
    for(auto edge:edges){
        g.addEdge(edge.first,edge.second);
    }
    return g.hasCycle();
}
/**
 * @brief static function to determine the shortest cost to traverse a graph
 *      using Dijkstra's algorithm
 * @param grid  2D grid of integers
 * @return  integer representing the shortest cost to traverse the graph
 */
int shortest_path(vector<vector<int> > grid){
    //return the shortest path len
    int n = grid.size();
    int m = grid[0].size();
    int size = n * m;
    Graph graph(grid, size, true);
    int total = graph.dijkstra(0, size - 1) + grid[0][0];
    return total;
}
/**
 * @brief static function to determine the biggest island in a graph. GIven a
 * two dimensional grid, containing only 0's and 1's. Each 1 represents land and
 * 0 represents water. The adjacent 1's form an island. Each land piece (x,y) is
 * connected to it's 4 neighbours (N,S,E,W). This function should return the
 * size of the largest island in the grid and 0 if there are no islands.
 * @param grid
 * @return
 */
