//
// Created by Ryan.Zurrin001 on 1/11/2022.
//

#ifndef GRAPH_CPP_UDEMY_H
#define GRAPH_CPP_UDEMY_H
#pragma ide diagnostic ignored "bugprone-branch-clone"
#pragma ide diagnostic ignored "readability-use-anyofallof"
#pragma ide diagnostic ignored "misc-no-recursion"
#include <bits/stdc++.h>
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
    [[nodiscard]] int degree() const {
        return neighbors.size();
    }
    bool hasEdge(GraphNode<T>* node1, GraphNode<T>* node2) {
        for (auto& neighbor : node1->neighbors) {
            if (neighbor.first == node2->data) {
                return true;
            }
        }
        return false;
    }
     void setParent(T parent_) { this->parent = parent_; }
     void setDistance(double distance_) { this->distance = distance_; }
     void setVisited(bool visited_) { this->visited = visited_; }

}; //___________________________end graph node class____________________________
/**
 * @brief a comparator for the priority queue that compares the weights of two
 *      nodes
 */
template<class T>
class  GraphNodeComparator {
public:
    // compare the weights of two nodes
    bool operator()(const pair<T, double> &a, const pair<T, double> &b) {
        return a.second > b.second;
    }
}; // ___________________end graph node comparator class________________________
// _______________________________begin graph class_____________________________
template <class T>
class  Graph {
     unordered_map<T, GraphNode<T>*> nodes;

     int V{};
     int E{};
     bool isDirected{};
     bool isWeighted{};
     void bfsUtil(T src,const bool* visited,const bool* recur);
     void dfsUtil( T src, bool* visited, bool* recur);
     void dfsUtil(T src, bool* visited);
     bool hasCycleUtil( T node, bool* visited, bool* recur,T parent);
     bool isConnectedUtil( int i, bool* visited);
     bool findComponentsUtil( int i, bool* visited, vector<T>& components);
     int traverse( T u, bool* visited);
     bool cycleFromVertexUtil(T node,  bool* visited);
     void findBackEdgesUtil( T node,bool* visited,bool* recur,vector<pair<T,T>>& backEdges);
     bool isBipartiteUtil( T node,  const bool* visited,int* color);
     int encode(int i, int j, int rows, int cols, bool rowMajor = true);
     double pathLengthUtil(T node1, T node2);

public:
     explicit Graph(bool isDirected_ = false);
     explicit Graph(vector<vector<T>> grid,int numbVertices,bool isDirected = false);
     explicit Graph(int vertices, bool isDirected = false);
     Graph(int vertices, pair<T,T> edges, bool isDirected = false);
     Graph(int vertices, tuple<T,T,double> w_edgs, bool isDirected = false);
     explicit Graph( vector<T> nodes, bool isDirected_ = false);
     explicit Graph(const string& fileName,char delimiter = ' ',
                    bool isDirected_ = false, bool isWeighted_ = false);
     Graph(T* nodes, int size, bool isDirected_ = false);
     void addVertex(T i);
     void addEdge(T node1, T node2);
     void addEdge(T node1, T node2, double weight);
     void removeEdge(T node1, T node2);
     void removeNode(T node);
     void bfs( T src);
     void dfs(T src);
     void findBackEdges(vector<pair<T,T>>& backEdges);
     void print();
     void printAllGraphData();
     bool isConnected();
     bool hasCycle();
     bool cycleFromVertex(T node);
     bool isBipartite();
     bool hasEdge(T node1, T node2);
     bool directed() const;
     int getV() const;
     int getE() const;
     double getWeight(T node1, T node2);
     int getNumberOfComponents();
     int degree(T node);
     T getId(T node);
     int minDistance(const int *pInt, const bool *pBoolean);
     int getCombinations();
     double pathLength(T node1, T node2);
     unordered_map<T, GraphNode<T>*> getNodes() const;
     vector<pair<T,T>> dijkstra(T src, T dest, bool print);
     T dijkstra(T src, T dest);
     vector<T> shortestPath(T src, T dest, bool print = false);
     vector<vector<T>> shortestPaths(T src, bool print = false);
     vector<vector<T>> getComponents();

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
    V = 0;
    E = 0;
    isWeighted = false;
} // ____________________________end graph constructor__________________________

template<class T>
Graph<T>::Graph(vector<vector<T>> grid, int numbVertices, bool isDirected) {
    // create graph with grid
    int rows = grid.size();
    int cols = grid[0].size();
    this->isDirected = isDirected;
    this->isWeighted = true;
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
                    if (!hasEdge(encodedIJ, encodedXY)) {
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
template<class T>
Graph<T>::Graph(int vertices, bool isDirected_) {
    isDirected = isDirected_;
    isWeighted = false;
    for (int i = 0; i < vertices; i++) {
        addVertex(i);
    }
} // ____________________________end graph constructor__________________________


template<class T>
Graph<T>::Graph(int vertices, pair<T, T> edges, bool isDirected) {
    isDirected = isDirected;
    isWeighted = false;
    for (int i = 0; i < vertices; i++) {
        addVertex(i);
    }
    for (int i = 0; i < edges.first.size(); i++) {
        addEdge(edges.first[i], edges.second[i]);
    }
} // ____________________________end graph constructor__________________________

template<class T>
Graph<T>::Graph(int vertices, tuple<T, T, double> w_edgs, bool isDirected) {
    isDirected = isDirected;
    isWeighted = true;
    for (int i = 0; i < vertices; i++) {
        addVertex(i);
    }
    for (int i = 0; i < get<0>(w_edgs).size(); i++) {
        addEdge(get<0>(w_edgs)[i], get<1>(w_edgs)[i], get<2>(w_edgs)[i]);
    }
} // ____________________________end graph constructor__________________________

/**
 * @brief Graph constructor with vector of nodes
 * @param nodes  : vector of nodes
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>
Graph<T>::Graph( vector<T> nodes,  bool isDirected_) {
    isDirected = isDirected_;
    isWeighted = false;
    for (int i = 0; i < nodes.size(); i++) {
        addVertex(nodes[i]);
    }
} // ____________________________end graph constructor__________________________
/**
 * @brief Graph constructor with a file name and delimiter
 * @param fileName  : file name
 * @param delimiter  : delimiter
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>
Graph<T>::Graph(const string& fileName,
                 char delimiter,
                bool isDirected_, bool isWeighted_) {
    isDirected = isDirected_;
    isWeighted = isWeighted_;
    ifstream file(fileName);
    if (!file.is_open()) {
        cout << "File not found" << endl;
        return;
    }
    // if weighted graph, read in the weights or else read in the edges
    if (isWeighted) {
        string line;
        while (getline(file, line)) {
            stringstream ss(line);
            string token;
            vector<string> tokens;
            while (getline(ss, token, delimiter)) {
                tokens.push_back(token);
            }
            addEdge(stoi(tokens[0]), stoi(tokens[1]), stod(tokens[2]));
        }
    } else {
        string line;
        while (getline(file, line)) {
            stringstream ss(line);
            string token;
            vector<string> tokens;
            while (getline(ss, token, delimiter)) {
                tokens.push_back(token);
            }
            addEdge(stoi(tokens[0]), stoi(tokens[1]));
        }
    }
} // ____________________________end graph constructor__________________________
/**
 * @brief Graph constructor with an array of nodes
 * @param nodes  : array of nodes
 * @param size  : size of array
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>
Graph<T>::Graph(T *nodes_, int size, bool isDirected_) {
    isDirected = isDirected_;
    isWeighted = false;
    for (int i = 0; i < size; i++) {
        addVertex(i);
        this->nodes[i]->data = nodes_[i];
    }
} // ____________________________end graph constructor__________________________

/**
 * @brief adds a vertex to the graph
 * @param i  : node to be added
 */
template<class T>
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
template<class T>
void Graph<T>::addEdge(T node1, T node2, double weight) {
    isWeighted = true;
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        throw invalid_argument("One of the nodes does not exist");
    }
    nodes[node1]->neighbors.push_back(make_pair(node2, weight));
    if (!isDirected) {
        nodes[node2]->neighbors.push_back(make_pair(node1, weight));
    }
    E++;
} // _____________________________end addEdge___________________________________
/**
 * @brief removes an edge from the graph
 * @param node1  : first node
 * @param node2  : second node
 */
template<class T>
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
template<class T>
void Graph<T>::removeNode( T node) {
    if (nodes.find(node) == nodes.end()) {
        throw invalid_argument("Node does not exist");
    }
    for ( auto neighbor : nodes[node]->neighbors) {
        removeEdge(node, neighbor.first);
    }
    nodes.erase(node);
    V--;
} // ________________________end removeNode_____________________________________
/**
 * @brief breadth first search
 * @param src  : source node to start search
 */
template<class T> 
void Graph<T>::bfs( T src) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    queue<T> q;
     bool* visited = new bool[V]{false};
    q.push(src);
    visited[src] = true;
    while (!q.empty()) {
        T node = q.front();
        q.pop();
        cout << node << " ";
        for ( auto neighbor : nodes[node]->neighbors) {
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
template<class T> 
void Graph<T>::dfs( T src) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
     bool* visited = new bool[V]{false};
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
template<class T> 
void Graph<T>::print() {
    for ( auto node : nodes) {
        cout << node.first << ": ";
        for ( auto neighbor : node.second->neighbors) {
            cout << neighbor.first << " ";
        }
        cout << endl;
    }
} // ________________________end print__________________________________________
/**
 * @brief prints the graph in adjacency matrix format with weights
 */
template<class T> 
void Graph<T>::printAllGraphData() {
    cout << "V: " << V << endl;
    cout << "E: " << E << endl;
    // print the weights of the edges as well
    for ( auto node : nodes) {
        cout << node.first << ":";
        for ( auto neighbor : node.second->neighbors) {
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
     bool* visited;
    T parent = -1;
    // set all nodes to unvisited
    visited = new bool[V];
     bool *recStack = new bool[V];
    for ( int i = 0; i < V; i++) {
        recStack[i] = false;
        visited[i] = false;
    }
    if (isDirected) {
        for ( auto node : nodes) {
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
 bool Graph<T>::cycleFromVertex(T node) {
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
    for ( auto node : nodes) {
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
bool Graph<T>::hasEdge(T node1, T node2) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        throw invalid_argument("Node does not exist");
    }
    // iterate over the neighbors of node1
    for ( auto neighbor : nodes[node1]->neighbors) {
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
 bool Graph<T>::directed() const {
    return this->isDirected;
} // ________________________end directed_______________________________________
/**
 * @brief number of vertices in the graph
 * @return  : number of vertices in the graph
 */
template<class T>
int Graph<T>::getV() const {
    return V;
} // ________________________end getV___________________________________________
/**
 * @brief number of edges in the graph
 * @return  : number of edges in the graph
 */
template<class T>
int Graph<T>::getE() const {
    return E;
} // ________________________end getE___________________________________________


/**
 * @brief gets the minimum distance between a node and all other nodes and
 *  returns the shortest path distance
 * @param pInt : pointer to an integer to store the shortest path distance
 * @param pBoolean : pointer to a boolean to store the result of the function
 * @return   the min distance between the node and all other nodes
 */
template<class T>
int Graph<T>::minDistance(const int *pInt, const bool *pBoolean) {
    int min = INT_MAX;
    int min_index;
    for (int v = 0; v < V; v++) {
        if (pBoolean[v] == false && pInt[v] <= min) {
            min = pInt[v];
            min_index = v;
        }
    }
    return min_index;
} // ________________________end minDistance____________________________________
/**
 * @brief counts the total combination pairs that can be made from the graph
 * @param C  : the size of the combinations to be counted
 * @return
 */
template<class T>
int Graph<T>::getCombinations() {
    bool *visited = new bool[V]{false};
    int total = V*(V - 1) / 2;
    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            auto combos = traverse(i, visited);
            total -= combos * (combos - 1) / 2;
        }
    }
    return total;
} // ________________________end getCombonations_______________________________

/**
 * @brief returns the length between vertices node1 and node2
 * @param node1  : first node
 * @param node2  : second node
 * @return  : length between the two nodes
 */
template<class T>
double Graph<T>::pathLength(T node1, T node2) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        throw invalid_argument("Node does not exist");
    }
    return pathLengthUtil(node1, node2);
} // ________________________end pathLength_____________________________________

/**
 * @brief gets the weight of an edge between two nodes
 * @param node1 : first node
 * @param node2 : second node
 * @return : weight of the edge
 */
template<class T>
double Graph<T>::getWeight( T node1,  T node2) {
    // nodes are out of bounds
    if (node1 >= V || node2 >= V) {
        return -1;
    }
    for ( auto neighbor : nodes[node1]->neighbors) {
        if (neighbor.first == node2) {
            return neighbor.second;
        }
    }
    return -1;
} // ________________________end getWeight_____________________________________

/**
 * @brief getter function for the nodes of the graph
 * @return  : unordered map of all the nodes in the graph
 */
template<class T>
 unordered_map<T, GraphNode<T> *> Graph<T>::getNodes() const {
    return nodes;
}// _________________________end getNodes_______________________________________

/**
 *  @brief algorithm to find the shortest path between two nodes
 * @param src  : source node
 * @param dest  : destination node
 * @return  : the shortest path between the two nodes
 */
template<class T>
T Graph<T>::dijkstra( T src,
                      T dest) {
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
         auto it = s.begin();
         auto curNode = it->second;
         auto distTilNow = it->first;
        s.erase(it);
        for ( auto node : nodes[curNode]->neighbors) {
             auto neighbor = node.first;
             auto currEdge = node.second;
            if (distTilNow + currEdge < distances[neighbor]) {
                 auto f = s.find({distances[neighbor], neighbor});
                if (f != s.end()) {
                    s.erase(f);
                }
                distances[neighbor] = distTilNow + currEdge;
                s.insert({distances[neighbor], neighbor});
            }
        }
    }
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
template<class T>
vector<pair<T,T>> Graph<T>::dijkstra( T src,
                                      T dest,
                                      bool print) {
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
         auto it = s.begin();
         auto curNode = it->second;
         auto distTilNow = it->first;
        s.erase(it);
        for ( auto node : nodes[curNode]->neighbors) {
             auto neighbor = node.first;
             auto currEdge = node.second;
            if (distTilNow + currEdge < distances[neighbor]) {
                 auto f = s.find({distances[neighbor], neighbor});
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
        for ( auto p : path) {
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
template<class T> 
vector<T> Graph<T>::shortestPath( T src,
                                  T dest,
                                  bool print) {
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
     bool *visited = new bool[V];
     bool *recStack = new bool[V];
    for ( int i = 0; i < V; i++) {
        visited[i] = false;
        recStack[i] = false;
    }
    queue<T> q;
     int* dist = new int[V]{0};
     T* prev = new T[V]{-1};
    q.push(src);
    visited[src] = true;
    dist[src] = 0;
    while (!q.empty()) {
        T u = q.front();
        q.pop();
        for ( auto neighbor : nodes[u]->neighbors) {
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
            for ( auto node : path) {
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
template<class T>
vector<vector<T>> Graph<T>::shortestPaths( T src,
                                           bool print) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    // determine if the graph is weighted
     bool weighted = false;
    for ( auto node : nodes) {
        for ( auto neighbor : node.second->neighbors) {
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
     bool *visited = new bool[V];
     bool *recStack = new bool[V];
    for ( int i = 0; i < V; i++) {
        visited[i] = false;
        recStack[i] = false;
    }
    queue<T> q;
     int* dist = new int[V]{0};
     T* prev = new T[V]{-1};

    q.push(src);
    visited[src] = true;
    dist[src] = 0;
    while (!q.empty()) {
        T u = q.front();
        q.pop();
        for ( auto neighbor : nodes[u]->neighbors) {
            if (!visited[neighbor.first]) {
                visited[neighbor.first] = true;
                dist[neighbor.first] = dist[u] + 1;
                prev[neighbor.first] = u;
                q.push(neighbor.first);
            }
        }
    }
    for ( int i = 0; i < V; i++) {
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
                for ( auto node : path) {
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

/**
 * @brief method to find all the connected nodes of each component in the graph
 *  and store them in a vector of vectors
 * @return vector of vectors of connected nodes
 */
template<class T>
vector<vector<T>> Graph<T>::getComponents() {
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++) {
        visited[i] = false;
    }
    vector<vector<T>> components;
    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            vector<T> component;
            findComponentsUtil(i, visited, component);
            components.push_back(component);
        }
    }
    delete[] visited;
    return components;
} // ________________________end findAllComponents______________________________


//##############################################################################
//***************************PRIVATE METHODS************************************
//##############################################################################

// ___________________________bfsUtil___________________________________________
/**
 * @brief private helper function to help with bfs traversals
 * @param src  source node
 * @param visited  vector of visited nodes
 * @param recur  vector representing a stack of nodes to keep track of what node
 *  a node was visited from
 */
template<class T> 
void Graph<T>::bfsUtil(
         T src,  const bool* visited,
         const bool* recur) {
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
        for ( auto neighbor : nodes[node]->neighbors) {
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
void Graph<T>::dfsUtil( T src,
                        bool *visited,
                        bool *recur) {
    visited[src] = true;
    recur[src] = true;
    for ( auto neighbor : nodes[src]->neighbors) {
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
template<class T> 
void Graph<T>::dfsUtil(  T src,
                        bool* visited) {
    visited[src] = true;
    for ( auto neighbor : nodes[src]->neighbors) {
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
template<class T>
bool Graph<T>::hasCycleUtil( T node,
                             bool* visited,
                             bool* recur,
                             T parent) {
    if (isDirected) {
        if (!visited[node]) {
            visited[node] = true;
            recur[node] = true;
            // recur all the vertices adjacent to this vertex
            for ( auto neighbor : nodes[node]->neighbors) {
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
        for ( auto neighbor : nodes[node]->neighbors) {
            if (!visited[neighbor.first]) {
                 bool cycle = hasCycleUtil(neighbor.first, visited, recur, node);
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
bool Graph<T>::isConnectedUtil( int i,
                                bool* visited) {
    visited[i] = true;
    if (!isDirected) { // undirected
        for ( auto neighbor : nodes[i]->neighbors) {
            if (!visited[neighbor.first]) {
                 bool cycle =
                        isConnectedUtil(neighbor.first, visited);
                if (cycle) {
                    return true;
                }
            }
        }
        return false;
    } else { // directed
        traverse(i, visited);
        for ( auto neighbor : nodes[i]->neighbors) {
            if (!visited[neighbor.first]) {
                return false;
            }
        }
        return true;
    }
} // ________________________end isConnectedUtil________________________________

/**
 * @brief private helper function to help in finding all the components of each
 *  connected component.
 * @param i node to start the search from
 * @param visited vector of visited nodes
 * @param components vector of components
 * @return true if the graph is connected, false otherwise
 */
template<class T>
bool
Graph<T>::findComponentsUtil(int i, bool *visited, vector<T> &components) {
    visited[i] = true;
    if (!isDirected) { // undirected
        for ( auto neighbor : nodes[i]->neighbors) {
            if (!visited[neighbor.first]) {
                findComponentsUtil(neighbor.first, visited, components);
            }
        }
        components.push_back({i});
        return false;
    } else { // directed
        traverse(i, visited);
        for ( auto neighbor : nodes[i]->neighbors) {
            if (!visited[neighbor.first]) {
                return false;
            }
        }
        components.push_back({i});
        return true;
    }
} // ________________________end isConnectedUtil________________________________


/**
 * @brief basic traversal of the graph
 * @param u node to look at in the visited vector
 * @param visited  vector of visited nodes
 */
template<class T>
int Graph<T>::traverse( T u,  bool *visited) {
    visited[u] = true;
    int count = 1;
    for ( auto neighbor : nodes[u]->neighbors) {
        if (!visited[neighbor.first]) {
            count += traverse(neighbor.first, visited);
        }
    }
    return count;
} // ________________________end traverse_______________________________________

/**
 * private utility function to determine if graph has a cycle from a given node
 * @param node  node to start the search from
 * @param visited  vector of visited nodes
 * @param parent  parent node of the current node
 */
template<class T>
bool Graph<T>::cycleFromVertexUtil( T node,
                                    bool *visited) {
    visited[node] = true;
    for ( auto neighbor : nodes[node]->neighbors) {
        if (!visited[neighbor.first]) {
            if (cycleFromVertexUtil(neighbor.first, visited)) {
                return true;
            }
        } else if (neighbor.first != node) {
            return true;
        }
    }
    return false;
}// ________________________end cycleFromVertexUtil____________________________
/**
 * @brief private utility function to determine if graph has a cycle
 * @param node  node to start the search from
 * @param visited  vector of visited nodes
 * @param parent  parent node of the current node
 */
template<class T>
void Graph<T>::findBackEdgesUtil( T node,
                                  bool *visited,
                                  bool *recur,
                                  vector<pair<T, T>> &backEdges) {
    visited[node] = true;
    recur[node] = true;
    for ( auto neighbor : nodes[node]->neighbors) {
        if (!visited[neighbor.first]) {
            findBackEdgesUtil(neighbor.first, visited, recur, backEdges);
        } else if (recur[neighbor.first]) {
            backEdges.push_back(make_pair(node, neighbor.first));
        }
    }
    recur[node] = false;
} // ________________________end findBackEdgesUtil_____________________________


/**
 * @brief if graph is disconnected will return the value of the total number of
 *  components in the graph
 * @return  number of components in the graph
 */
template<class T>
int Graph<T>::getNumberOfComponents() {
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++) {
        visited[i] = false;
    }
    int count = 0;
    for (int i = 0; i < V; i++) {
        if (!visited[i]) {
            count++;
            isConnectedUtil(i, visited);
        }
    }
    delete[] visited;
    return count;
} // ________________________end getNumberOfComponents__________________________

template<class T>
int Graph<T>::degree(T node) {
    return nodes[node]->neighbors.size();
}


/**
 * @brief method to get the id of a graph node
 * @param node : node to get the id of
 * @return  id of the node or -1 if node is not in the graph
 */
template<class T>
T Graph<T>::getId(T node) {
    auto it = nodes.find(node);
    if (it != nodes.end()) {
        return it->first;
    }
    return -1;
} // ________________________end getId_________________________________________
/**
 * @brief method to determine if the graph is bipartite
 * @param node : node to start the search from
 * @param visited : vector of visited nodes
 * @param color  : color array to store the color of each node
 * @return  true if the graph is bipartite, false otherwise
 */
template<class T>
bool Graph<T>::isBipartiteUtil( T node,
                                const bool *visited,
                                int *color) {
    color[node] = 1;
    visited[node] = true;
    queue<T> q;
    q.push(node);
    while (!q.empty()) {
        auto u = q.front();
        q.pop();
        for ( auto neighbor : nodes[u]->neighbors) {
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
} // ________________________end isBipartiteUtil________________________________
/**
 * @brief private utility function to decode the i and j values to either
 * row or column major order index
 * @param i : row index
 * @param j  : column index
 * @param rows  : number of rows
 * @param cols  : number of columns
 * @param rowMajor  : true if row major order is used else false
 * @return  single index of the ith row and jth column
 */
template<class T>
int Graph<T>::encode(int i, int j, int rows, int cols, bool rowMajor) {
    if (rowMajor) {
        return i * cols + j;
    } else {
        return i + j * rows;
    }
}

template<class T>
double Graph<T>::pathLengthUtil(T node1, T node2) {
    if (node1 == node2) {
        return 0;
    }
    bool *visited = new bool[V];
    for (int i = 0; i < V; i++) {
        visited[i] = false;
    }
    auto *dist = new double[V];
    for (int i = 0; i < V; i++) {
        dist[i] = INT_MAX;
    }
    if (isWeighted) {
        dist[node1] = 0;
        queue<T> q;
        q.push(node1);
        while (!q.empty()) {
            auto u = q.front();
            q.pop();
            for ( auto neighbor : nodes[u]->neighbors) {
                if (!visited[neighbor.first] && dist[neighbor.first] > dist[u] + neighbor.second) {
                    dist[neighbor.first] = dist[u] + neighbor.second;
                    q.push(neighbor.first);
                }
            }
            visited[u] = true;
        }
    } else {
        dist[node1] = 0;
        queue<T> q;
        q.push(node1);
        while (!q.empty()) {
            auto u = q.front();
            q.pop();
            for ( auto neighbor : nodes[u]->neighbors) {
                if (!visited[neighbor.first] && dist[neighbor.first] > dist[u] + 1) {
                    dist[neighbor.first] = dist[u] + 1;
                    q.push(neighbor.first);
                }
            }
            visited[u] = true;
        }
    }
    delete[] visited;
    return dist[node2];
}


void dfs_utl(vector<vector<int>> mat, vector<vector<bool>>& visited,
             vector<vector<int>>& mem, int i, int j, int m, int n) {

    visited[i][j] = 1;

    int dx[] = {-1,1,0,0};
    int dy[] = {0,0,1,-1};

    int cnt = 0;
    for(int k=0;k<4;k++){
        int nx = i + dx[k];
        int ny = j + dy[k];

        if(nx>=0 and ny>=0 and nx<m and ny<n and mat[nx][ny]>mat[i][j]){
            int subProblemCnt = 0;
            if(visited[nx][ny]){
                cnt = max(cnt,1+mem[nx][ny]);
            }
            else{
                dfs_utl(mat,visited,mem,nx,ny,m,n);
                cnt = max(cnt,1+mem[nx][ny]);
            }
        }
    }
    mem[i][j] = cnt;
}

int longestPathSequence(vector<vector<int> > matrix){
    //return the longest path sequence
    int m = matrix.size();
    int n = matrix[0].size();
    vector<vector<bool> > visited(m+1,vector<bool>(n+1,0));
    vector<vector<int> > cache(m+1,vector<int>(n+1,0));
    int ans = 0;
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            dfs_utl(matrix,visited,cache,i,j,m,n);
            ans = max(ans,cache[i][j]);
        }
    }
    return ans+1;
}


int count_pairs(int N, const vector<pair<int,int> >& astronauts){
    //return the number of pairs
    Graph<int> g;
    for (int i = 0; i < N; i++) {
        g.addVertex(i);
    }
    for(auto edge:astronauts){
        g.addEdge(edge.first,edge.second);
    }
    return g.getCombinations();
}
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
static int shortest_path(vector<vector<int> > grid){
    //return the shortest path len
    auto n = grid.size();
    auto m = grid[0].size();
    int size = n * m;
    Graph graph(grid, size, true);
    int total = graph.dijkstra(0, size - 1) + grid[0][0];
    return total;
}


#endif //GRAPH_CPP_UDEMY_H
