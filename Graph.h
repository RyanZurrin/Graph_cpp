#pragma ide diagnostic ignored "bugprone-branch-clone"
#pragma ide diagnostic ignored "readability-use-anyofallof"
#pragma ide diagnostic ignored "misc-no-recursion"
#pragma ide diagnostic ignored "hicpp-use-auto"
#ifndef GRAPH_CPP_GRAPH_H
#define GRAPH_CPP_GRAPH_H
#include <bits/stdc++.h>
using namespace std;

template<typename GraphNode>
class GraphNodeIterator {
public:
    using ValueType = typename GraphNode::ValueType;
    using PointerType = ValueType*;
    using ReferenceType = ValueType&;
    using ConstPointerType = const ValueType*;
    using ConstReferenceType = const ValueType&;
public:
    explicit GraphNodeIterator(GraphNode* node) : node_(node) {}
    GraphNodeIterator(const GraphNodeIterator& other) : node_(other.node_) {}
    GraphNodeIterator& operator=(const GraphNodeIterator& other) {
        if (this != &other) {
            node_ = other.node_;
        }
        return *this;
    }
    bool operator==(const GraphNodeIterator& other) const {
        return node_ == other.node_;
    }
    bool operator!=(const GraphNodeIterator& other) const {
        return node_ != other.node_;
    }
    GraphNodeIterator& operator++() {
        node_ = node_->next_;
        return *this;
    }
    GraphNodeIterator operator++(int) {
        GraphNodeIterator temp(*this);
        node_ = node_->next_;
        return temp;
    }
    ReferenceType operator*() {
        return node_->value_;
    }
    ConstPointerType operator->() const {
        return &(node_->value_);
    }
private:
    GraphNode* node_;

};
template<typename Graph>
class GraphIterator {
public:
    using ValueType = typename Graph::ValueType;
    using PointerType = ValueType*;
    using ReferenceType = ValueType&;
    using GraphNodePtr = typename Graph::GraphNodePtr;
    using ConstPointerType = const ValueType*;
    using ConstReferenceType = const ValueType&;
public:
    explicit GraphIterator(GraphNodePtr node) : node_(node) {}
    GraphIterator(const GraphIterator& other) : node_(other.node_) {}
    GraphIterator& operator=(const GraphIterator& other) {
        node_ = other.node_;
        return *this;
    }
    bool operator==(const GraphIterator& other) const {
        return node_ == other.node_;
    }
    bool operator!=(const GraphIterator& other) const {
        return node_ != other.node_;
    }
    bool operator>(const GraphIterator& other) const {
        return node_ > other.node_;
    }
    bool operator<(const GraphIterator& other) const {
        return node_ < other.node_;
    }
    GraphIterator& operator++() {
        node_ = node_->next_;
        return *this;
    }
    GraphIterator operator++(int) {
        GraphIterator temp(*this);
        node_ = node_->next_;
        return temp;
    }
    ReferenceType operator*() {
        return node_->value_;
    }
    ConstPointerType operator->() const {
        return &(node_->value_);
    }
private:
    GraphNodePtr node_;
};
// ___________________________begin graph node class____________________________
/**
 * @brief The GraphNode class
 * @tparam T  the type of the node
 */
template <class T>
class GraphNode {
public:
    using ValueType = T;
    using Iterator = typename vector<GraphNode<T>*>::iterator;
    using ConstIterator = typename vector<GraphNode<T>*>::const_iterator;
    using ReverseIterator = typename vector<GraphNode<T>*>::reverse_iterator;
    using ConstReverseIterator = typename vector<GraphNode<T>*>::const_reverse_iterator;
    using GraphNodePtr = GraphNode<T>*;
    using GraphNodePtrConst = const GraphNode<T>*;
    using GraphNodePtrConstPtr = const GraphNode<T>**;
public:
    T data;
    list<pair<int, double>> neighbors;
    std::shared_ptr<T> parent;
    bool visited{};
    double distance{};
    GraphNode() : data(0), parent(nullptr), visited(false), distance(0) {}
    explicit GraphNode(T data) {
        this->data = data;
        this->visited = false;
        this->distance = 0;
        this->parent = std::make_shared<T>(0);
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
    // get edge
    double getEdge(GraphNode<T>* node1, GraphNode<T>* node2) {
        for (auto& neighbor : node1->neighbors) {
            if (neighbor.first == node2->data) {
                return neighbor.second;
            }
        }
        return 0;
    }
     void setParent(T parent_) { this->parent = parent_; }
     void setDistance(double distance_) { this->distance = distance_; }
     void setVisited(bool visited_) { this->visited = visited_; }
    // ________________________Iterator Methods_________________________________
    Iterator begin() { return neighbors.begin(); }
    Iterator end() { return neighbors.end(); }
    [[nodiscard]] ConstIterator begin() const { return neighbors.begin(); }
    [[nodiscard]] ConstIterator end() const { return neighbors.end(); }
    [[nodiscard]] ConstIterator cbegin() const { return neighbors.cbegin(); }
    [[nodiscard]] ConstIterator cend() const { return neighbors.cend(); }
    ReverseIterator rbegin() { return neighbors.rbegin(); }
    ReverseIterator rend() { return neighbors.rend(); }
    [[nodiscard]] ConstIterator crbegin() const { return neighbors.rbegin(); }
    [[nodiscard]] ConstIterator crend() const { return neighbors.rend(); }
    // ______________________End Iterator Methods_______________________________
    ~GraphNode() {
        neighbors.clear();
    }
    const auto &toString() {
        // find node in list and print the data of the node
        for (auto& neighbor : neighbors) {
            if (neighbor.first == data) {
                return neighbor.first;
            }
        }
    }
    const auto &getData() {
        return data;
    }
};
//______________________________end graph node class____________________________
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
};
// ___________________end graph node comparator class___________________________
// _______________________________begin graph class_____________________________
template <class T>
class  Graph {
public:
    // __________________________Graph Constructors_____________________________
     explicit Graph(bool isDirected_ = false);

     explicit Graph(vector<vector<T>> grid, int numbVertices, bool isDirected = false);

     explicit Graph(int vertices, bool isDirected = false);

     Graph(int vertices, pair<T,T> edges, bool isDirected = false);

     Graph(int vertices, tuple<T,T,double> w_edgs, bool isDirected = false);

     explicit Graph( vector<T> nodes, bool isDirected_ = false);

     explicit Graph(const string& fileName, char delimiter = ' ',
                    bool isDirected_ = false,
                    bool isWeighted_ = false);

     Graph(T* nodes, int size, bool isDirected_ = false);

     Graph(const Graph& other);

     Graph(Graph&& other);

     Graph& operator=(const Graph& other);

     Graph& operator=(Graph&& other);

     // __________________________Graph Destructor______________________________
     ~Graph();

     // __________________________Graph Setter Methods__________________________
     void addVertex(T i); // add a vertex to the graph
     void addVertex(pair<int,T> i); // add a vertex to the graph
     // add node to graph
     void addNode(GraphNode<T>* node);
     void addEdge(T node1, T node2); // add an edge to the graph
     void addEdge(T node1, T node2, double weight); // add weighted edge
     void addEdge(pair<int,T> node1, pair<int, T> node2); // add an edge to the graph
     void addEdge(pair<int,T> node1, pair<int, T> node2, double weight);
     void setData(T data_, int node); // set the data of a node
     void setWeight(int node1, int node2, double weight_); // set the weight of an edge

     // __________________________Graph Getter Methods__________________________
     std::shared_ptr<GraphNode<T>> getVertex(T node) const;
     [[nodiscard]] unordered_map<T, GraphNode<T>*> getNodes() const;
     [[nodiscard]] int  getV() const;
     [[nodiscard]] int  getE() const;
     vector<T>          getNeighbors(T node);
     vector<vector<T>>  getComponents();
     T                  getId(T node);
     double             getAvgPathLength();
     double             getAvgDegree();
     double             getAvgClusteringCoefficient();
     double             getWeight(T node1, T node2);
     int                getNumberOfComponents();
     int                getDegree(T node);
     int                getCombinations();
     vector<vector<T>>  getAdjMatrix();
     vector<pair<pair<T,T>,double>>  getEdges();
     vector<pair<T,T>>  getBridges();
     vector<T>          getArticulationPoints();
     Graph<T>           getTranspose();

     // ____________________Graph Removal Methods_______________________________
     void               removeEdge(T node1, T node2);
     void               removeVertex(T node);
     void               emptyGraph();

     // _____________________Graph Boolean Methods______________________________
     [[nodiscard]] bool directed() const;
     bool               connected();
     bool               cycle();
     bool               cycleFromVertex(T node);
     bool               bipartite();
     bool               hasEdge(T node1, T node2);
     bool               isTree();
     bool               isDAG();

     // ____________________Graph Traversal Algorithms_____________________________
     void               bfs( T src);
     void               dfs(T src);
     double             dijkstra(T src, T dest);
     vector<T>          shortestPath(T src, T dest, bool print = false);
     vector<vector<T>>  shortestPaths(T src, bool print = false);
     vector<pair<T,T>>  dijkstra(T src, T dest, bool print);
     // bellman ford
     vector<pair<T,T>>  bellmanFord(T src,bool print = false);

     // ______________________Graph Misc. Methods_______________________________
     void findBackEdges(vector<pair<T,T>>& backEdges);
     int minDistance(const int *pInt, const bool *pBoolean);
     double pathLength(T node1, T node2);
     vector<T> topologicalSort();
     vector<vector<T>> stronglyConnectedComponents();


     // ______________________Graph Print Methods_______________________________
     void print();
     void printAllGraphData();
     void printNodeData();
     void iteratorPrinter();

    //_________________________Iterator Class Methods____________________________
    using ValueType = T;
    using GraphIterator = typename unordered_map<T, GraphNode<T>*>::iterator;
    using ConstGraphIterator = typename unordered_map<T, GraphNode<T>*>::const_iterator;
    using GraphNodeIterator = typename GraphNode<T>::Iterator;
    using ConstGraphNodeIterator = typename GraphNode<T>::ConstIterator;
    using ReverseGraphNodeIterator = typename GraphNode<T>::ReverseIterator;
    using ConstReverseGraphNodeIterator = typename GraphNode<T>::ConstReverseIterator;
    using GraphNodeConstIterator = typename GraphNode<T>::ConstIterator;
    GraphIterator begin() { return nodes.begin(); }
    GraphIterator end() { return nodes.end(); }
    ConstGraphIterator cbegin() const { return nodes.cbegin(); }
    ConstGraphIterator cend() const { return nodes.cend(); }
    GraphIterator find(T key) { return nodes.find(key); }
    ConstGraphIterator find(T key) const { return nodes.find(key); }



private:
    unordered_map<int, std::shared_ptr<GraphNode<T>>> nodes; // the graph
    int     V{}; // number of vertices
    int     E{}; // number of edges
    bool    isDirected{}; // is the graph directed
    bool    isWeighted{}; // is the graph weighted
    double  avgDegree{}; // average degree of the graph
    double  avgPathLenth{}; // average path length of the graph
    double  globalClusteringCoef{}; // global clustering coefficient of the graph
    int     deletedVertices{}; // number of vertices deleted from the graph
    // __________________________Graph Helper Methods___________________________
    double  calculateAverageDegree();
    double  calculateAveragePathLength();
    double  clusteringCoef(std::shared_ptr<GraphNode<T>>  node);
    double  calculateGlobalClusteringCoefficient();
    void    bfsUtil(T src, const bool* visited,const bool* recur);
    void    dfsUtil( T src, bool* visited, bool* recur);
    void    dfsUtil(T src, bool* visited);
    void    dfsUtil(int v, bool visited[], vector<T> &comp);
    bool    hasCycleUtil( T node, bool* visited, bool* recur, T parent);
    bool    isConnectedUtil( int i, bool* visited);
    bool    findComponentsUtil( int i, bool* visited, vector<T>& components);
    int     traverse( T u, bool* visited);
    bool    cycleFromVertexUtil(T node,  bool* visited);
    void    findBackEdgesUtil( T node, bool* visited, bool* recur,
                            vector<pair<T,T>>& backEdges);
    bool    isBipartiteUtil( T node,  bool* visited, int* color);
    int     encode(int i, int j, int rows, int cols, bool rowMajor = true);
    double  pathLengthUtil(T node1, T node2);
    void    fillOrder(int v, bool visited[], stack<T> &Stack);
    void    sccUtil(int u, int* disc, int* low, stack<T> *st, bool* stackMember,
                    vector<vector<T>> &scc);
    void    bridgeUtil(int u, bool* visited, int* disc, int* low, int* parent,
                    vector<pair<T,T>> &bridges);
    void    articulationPointUtil(int i, bool* visited, int* disc, int* low, int* parent,
                                  vector<T>& articulationPoints);



};
// _______________________________end graph class_______________________________
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
    avgDegree = 0;
    avgPathLenth = 0;
    globalClusteringCoef = 0;
    deletedVertices = 0;
}
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
    deletedVertices = 0;
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
    deletedVertices = 0;
}
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
    deletedVertices = 0;
}
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
    deletedVertices = 0;
}
/**
 * @brief Graph constructor with vector of nodes
 * @param nodes  : vector of nodes
 * @param isDirected_  : default is false, set to true if graph is directed
 */
template<class T>
Graph<T>::Graph( vector<T> nodes,
                 bool isDirected_) {
    isDirected = isDirected_;
    isWeighted = false;
    for (int i = 0; i < nodes.size(); i++) {
        addVertex(nodes[i]);
    }
}
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
    deletedVertices = 0;
}
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
    deletedVertices = 0;
}
/**
 * @brief copy constructor
 * @param other  : graph to be copied
 */
template<class T>
Graph<T>::Graph(const Graph &other) { // copy constructor
    isDirected = other.isDirected;
    isWeighted = other.isWeighted;
    deletedVertices = other.deletedVertices;
    for (int i = 0; i < other.nodes.size(); i++) {
        addVertex(other.nodes[i]->data);
    }
    // copy edges from other graph and if graph is weighted, copy weights
    for (int i = 0; i < other.nodes.size(); i++) {
        for (int j = 0; j < other.nodes[i]->edges.size(); j++) {
            if (isWeighted) {
                addEdge(i, other.nodes[i]->edges[j]->data,
                        other.nodes[i]->edges[j]->weight);
            } else {
                addEdge(i, other.nodes[i]->edges[j]->data);
            }
        }
    }
}
/**
 * @brief move constructor
 * @param other  : graph to be moved
 */
template<class T>
Graph<T>::Graph(Graph &&other) { // move constructor
    isDirected = other.isDirected;
    isWeighted = other.isWeighted;
    deletedVertices = other.deletedVertices;
    nodes = other.nodes;
    other.nodes.clear();
}
/**
 * @brief copy assignment operator
 * @param other  : graph to be copied
 * @return  : reference to this graph
 */
template<class T>
Graph<T> &Graph<T>::operator=(const Graph &other) { // copy assignment
    if (this != &other) {
        isDirected = other.isDirected;
        isWeighted = other.isWeighted;
        deletedVertices = other.deletedVertices;
        nodes.clear();
        for (int i = 0; i < other.nodes.size(); i++) {
            addVertex(other.nodes[i]->data);
        }
        // copy edges from other graph and if graph is weighted, copy weights
        for (int i = 0; i < other.nodes.size(); i++) {
            for (int j = 0; j < other.nodes[i]->edges.size(); j++) {
                if (isWeighted) {
                    addEdge(i, other.nodes[i]->edges[j]->data,
                            other.nodes[i]->edges[j]->weight);
                } else {
                    addEdge(i, other.nodes[i]->edges[j]->data);
                }
            }
        }
    }
    return *this;
}
/**
 * @brief move assignment operator
 * @tparam T
 * @param other  : graph to be moved
 * @return  : reference to this graph
 */
template<class T>
Graph<T> &Graph<T>::operator=(Graph &&other) { // move assignment
    if (this != &other) {
        isDirected = other.isDirected;
        isWeighted = other.isWeighted;
        deletedVertices = other.deletedVertices;
        nodes = other.nodes;
        other.nodes.clear();
    }
    return *this;
}

template<class T>
void Graph<T>::addVertex(pair<int, T> i) {
    // add a vertex to the graph
    if (nodes.find(i.first) != nodes.end()) {
        return;
    } else {
        nodes[i.first] = std::make_unique<GraphNode<T>>(i.first);
        nodes[i.first]->data = i.second;
    }
    V++;
}
/**
 * @brief adds a vertex to the graph
 * @param i  : node to be added
 */
template<class T>
void Graph<T>::addVertex(T i) {
    // check to see if the node already exists
    if (nodes.find(i) != nodes.end()) {
        return;
    } else {
        nodes[i] = std::make_shared<GraphNode<T>>(i);
        //nodes[i] = new GraphNode<T>(i);
        V++;
    }
}
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
    // check to see if the edge already exists
    for (auto &neighbor : nodes[node1]->neighbors) {
        if (neighbor.first == node2) {
            return;
        }
    }
    nodes[node1]->neighbors.push_back(make_pair(node2, 0));
    if (!isDirected) {
        nodes[node2]->neighbors.push_back(make_pair(node1, 0));
    }
    E++;
}
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
    // check to see if the edge already exists
    for (auto &neighbor : nodes[node1]->neighbors) {
        if (neighbor.first == node2) {
            return;
        }
    }
    nodes[node1]->neighbors.push_back(make_pair(node2, weight));
    if (!isDirected) {
        nodes[node2]->neighbors.push_back(make_pair(node1, weight));
    }
    E++;
}

template<class T>
void Graph<T>::addEdge(pair<int, T> node1, pair<int, T> node2) {
    if (nodes.find(node1.first) == nodes.end() || nodes.find(node2.first) == nodes.end()) {
        throw invalid_argument("One of the nodes does not exist");
    }
    // check to see if the edge already exists
    for (auto &neighbor : nodes[node1.first]->neighbors) {
        if (neighbor.first == node2.first) {
            return;
        }
    }
    nodes[node1.first]->neighbors.push_back(make_pair(node2.first, node2.second));
    if (!isDirected) {
        nodes[node2.first]->neighbors.push_back(make_pair(node1.first, node1.second));
    }
    E++;
}
template<class T>
void Graph<T>::addEdge(pair<int, T> node1, pair<int, T> node2, double weight) {
    isWeighted = true;
    if (nodes.find(node1.first) == nodes.end() || nodes.find(node2.first) == nodes.end()) {
        throw invalid_argument("One of the nodes does not exist");
    }
    // check to see if the edge already exists
    for (auto &neighbor : nodes[node1.first]->neighbors) {
        // if it is a string type then we need to compare the strings
        if (neighbor.first == node2.first) {
            return;
        }
    }
    nodes[node1.first]->neighbors.push_back(make_pair(node2.first, weight));
    if (!isDirected) {
        nodes[node2.first]->neighbors.push_back(make_pair(node1.first, weight));
    }
    E++;
}
template<class T>
void Graph<T>::setData(T data_, int node) {
    //check that the node exists
    if (nodes.find(node) == nodes.end()) {
        cout << "Node " << node << " does not exist" << endl;
        return;
    }
    nodes[node]->data = data_;
}

template<class T>
void Graph<T>::setWeight(int node1, int node2, double weight_) {
    //check that the node exists
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        cout << "Node " << node1 << " or " << node2 << " does not exist" << endl;
        return;
    }
    //check that the edge exists
    for (auto &neighbor : nodes[node1]->neighbors) {
        if (neighbor.first == node2) {
            cout << "Edge " << node1 << "-->" << node2 << " exist" << endl;
            cout << "setting weight from " << neighbor.second << " to " << weight_ << endl;
            neighbor.second = weight_;
            return;
        }
    }
    cout << "Edge " << node1 << "-->" << node2 << " does not exist" << endl;
}

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
}
/**
 * @brief removes a vertex from the graph
 * @param i  : node to be removed
 */
template<class T>
void Graph<T>::removeVertex( T node) {
    // check to see if the node exists
    if (nodes.find(node) == nodes.end()) {
        throw invalid_argument("The node does not exist");
    }
    // remove all edges from the node
    for (auto &neighbor : nodes[node]->neighbors) {
        cout << "removing edge" << neighbor.first << endl;
        nodes[neighbor.first]->neighbors.remove_if([&](pair<T, double> neighbor) {
            return neighbor.first == node;
        });
    }
    //delete nodes[node];
    nodes.erase(node);
    deletedVertices++;
    V--;
}

template<class T>
void Graph<T>::emptyGraph() {
    nodes.clear();
    V = 0;
    E = 0;
    deletedVertices = 0;
}

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
}
// ________________________________end bfs______________________________________
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
}
// _________________________________end dfs_____________________________________
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
// _________________________________end findBackEdges___________________________
/**
 * @brief prints the graph in adjacency list format
 */
template<class T> 
void Graph<T>::print() {
   // print using iterator
    for (auto node : nodes) {
        cout << node.first << ": ";
        for (auto neighbor : node.second->neighbors) {
            cout << neighbor.first << " ";
        }
        cout << endl;
    }
}
// ________________________end print____________________________________________
/**
 * @brief prints the graph in adjacency matrix format with weights
 */
template<class T> 
void Graph<T>::printAllGraphData() {
    avgDegree = calculateAverageDegree();
    avgPathLenth = calculateAveragePathLength();
    globalClusteringCoef = calculateGlobalClusteringCoefficient();
    cout << "V: " << V << endl;
    cout << "E: " << E << endl;
    // print the weights of the edges as well
    for ( auto node : nodes) {
        cout << node.first << ":";
        for ( auto neighbor : node.second->neighbors) {
            cout << "(" << neighbor.second << ")" << neighbor.first << ", ";
        }
        cout << endl;
    }
    // print out the average degree
    cout << "Average degree: " << avgDegree << endl;
    // print out the average path length
    cout << "Average path length: " << avgPathLenth << endl;
    // print out the global clustering coefficient
    cout << "Global clustering coefficient: " << globalClusteringCoef << endl;
}

template<class T>
void Graph<T>::printNodeData() {
    // use basic for loop to print out the data
    for (int i = 0; i < V; i++) {
        // check that the node exists
        if (nodes.find(i) != nodes.end()) {
            cout << "Node: " << i << endl;
            cout << "Degree: " << nodes[i]->neighbors.size() << endl;
            // print out the data of the node
            cout << "Data: " << nodes[i]->data << endl;
            cout << "Neighbors: ";
            for (auto neighbor : nodes[i]->neighbors) {
                cout << neighbor.first << " ";
            }
            cout << endl;
        }
    }
}

/**
 * @brief determines if the graph is connected
 * @return  : true if the graph is connected, false otherwise
 */
template<class T>
bool Graph<T>::connected() {
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
}
// _____________________________end isConnected_________________________________
/**
 * @brief determines if the graph has a cycle
 * @warning : this function does not work for disconnected undirected graphs
 * @return  : true if the graph has a cycle, false otherwise
 */
template<class T>
bool Graph<T>::cycle() {

    T parent = -1;
    // set all nodes to unvisited
    bool *visited = new bool[V];
     bool *recStack = new bool[V];
    for ( int i = 0; i < V; i++) {
        recStack[i] = false;
        visited[i] = false;
    }
    if (isDirected) {
        for ( const auto& node : nodes) {
            if (hasCycleUtil(node.first, visited, recStack, parent)) {
                delete[] recStack;
                delete[] visited;
                return true;
            }
        }
    } else {
        auto result = hasCycleUtil(0, visited, recStack, parent);
        delete[] recStack;
        delete[] visited;
        return result;
    }
    delete[] recStack;
    delete[] visited;
    return false;
}
// _____________________________end hasCycle____________________________________
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
}
// _____________________________end cycleFromVertex_____________________________
/**
 * @brief determines if the graph is bipartite
 * @return  : true if the graph is bipartite, false otherwise
 */
template<class T>
bool Graph<T>::bipartite() {
    int color[V];
    bool *visited = new bool[V];
    if (!isDirected) {
        for (int i = 0; i < V; i++) {
            visited[i] = false;
            color[i] = -1;
        }
        for ( auto node : nodes) {
            if (color[node.first] == -1) {
                if (!isBipartiteUtil(node.first, visited, color)) {
                    delete[] visited;
                    return false;
                }
            }
        }
    } else {
        for (int i = 0; i < V; i++) {
            visited[i] = false;
            color[i] = -1;
        }
        for ( auto node : nodes) {
            if (color[node.first] == -1) {
                if (!isBipartiteUtil(node.first, visited, color)) {
                    delete[] visited;
                    return false;
                }
            }
        }
    }
    delete[] visited;
    return true;
}
// _____________________________end isBipartite_________________________________
/**
 * @brief determines if an edge exists between two nodes
 * @param node1  : first node
 * @param node2  : second node
 * @return  : true if an edge exists, false otherwise
 */
template<class T>
bool Graph<T>::hasEdge(T node1, T node2) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        return false;
    }
    // iterate over the neighbors of node1
    for ( auto neighbor : nodes[node1]->neighbors) {
        if (neighbor.first == node2) {
            return true;
        }
    }
    return false;
}


template<class T>
bool Graph<T>::isTree() {
    // check if the graph is directed and has no cycles
    if (isDirected) {
        return false;
    }
    // check if the graph has cycles
    bool cycles = this->cycle();
    if (cycles) {
        return false;
    }
    return true;
}

template<class T>
bool Graph<T>::isDAG() {
    // check if the graph is directed and acyclic
    if (!isDirected) {
        return false;
    }
    // check if the graph has cycles
    bool cycles = this->cycle();
    if (cycles) {
        return false;
    }
    return true;
}

// _____________________________end findEdge____________________________________
/**
 * @brief returns the isDirected property of the graph
 * @return  : true if the graph is directed, false otherwise
 */
template<class T>
 bool Graph<T>::directed() const {
    return this->isDirected;
}
// ________________________end directed_________________________________________
/**
 * @brief number of vertices in the graph
 * @return  : number of vertices in the graph
 */
template<class T>
int Graph<T>::getV() const {
    return V;
}
// ________________________end getV_____________________________________________
/**
 * @brief number of edges in the graph
 * @return  : number of edges in the graph
 */
template<class T>
int Graph<T>::getE() const {
    return E;
}
template<class T>
std::shared_ptr<GraphNode<T>> Graph<T>::getVertex(T node) const {
    return nodes.at(node);
}
// __________________________end getVertex______________________________________
template<class T>
vector<T> Graph<T>::getNeighbors(T node) {
    vector<T> neighbors;
    // check that the node exists
    if (nodes.find(node) == nodes.end()) {
        cout << "Node "<< node << " does not exist" << endl;
        return neighbors;
    }
    for (auto it = nodes.at(node)->neighbors.begin(); it != nodes.at(node)->neighbors.end(); it++) {
        neighbors.push_back(it->first);
    }
    return neighbors;
}
// __________________________end getNeighbors___________________________________
// ________________________end getE_____________________________________________
template<class T>
double Graph<T>::getAvgPathLength() {
    avgPathLenth = calculateAveragePathLength();
    return avgPathLenth;
}
// ________________________end getAvgPathLength_________________________________
template<class T>
double Graph<T>::getAvgDegree() {
    avgDegree = calculateAverageDegree();
    return avgDegree;
}
// ________________________end getAvgDegree_____________________________________
template<class T>
double Graph<T>::getAvgClusteringCoefficient() {
    globalClusteringCoef = calculateGlobalClusteringCoefficient();
    return globalClusteringCoef;
}
// ________________________end getAvgClusteringCoefficient______________________
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
}
// ________________________end minDistance______________________________________
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
}
// ________________________end getCombinations_________________________________
template<class T>
vector<pair<pair<T,T>, double>> Graph<T>::getEdges() {
    // return the edges of the graph incuding weights between each node
    vector<pair<pair<T,T>, double>> edges;
    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        for (auto it2 = it->second->neighbors.begin(); it2 != it->second->neighbors.end(); it2++) {
            edges.push_back(make_pair(make_pair(it->first, it2->first), it2->second));
        }
    }
    return edges;
}
/**
 * @brief returns the length between vertices node1 and node2
 * @param node1  : first node
 * @param node2  : second node
 * @return  : length between the two nodes
 */
template<class T>
double Graph<T>::pathLength(T node1, T node2) {
    if (nodes.find(node1) == nodes.end() || nodes.find(node2) == nodes.end()) {
        return -1;
    }
    return pathLengthUtil(node1, node2);
}
// ________________________end pathLength_______________________________________

/**
 * @brief gets the weight of an edge between two nodes
 * @param node1 : first node
 * @param node2 : second node
 * @return : weight of the edge
 */
template<class T>
double Graph<T>::getWeight( T node1,
                            T node2) {
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
}
// ________________________end getWeight________________________________________

/**
 * @brief getter function for the nodes of the graph
 * @return  : unordered map of all the nodes in the graph
 */
template<class T>
unordered_map<T, GraphNode<T> *> Graph<T>::getNodes() const {
    return nodes;
}
// _________________________end getNodes________________________________________

/**
 *  @brief algorithm to find the shortest path between two nodes
 * @param src  : source node
 * @param dest  : destination node
 * @return  : the shortest path between the two nodes
 */
template<class T>
double Graph<T>::dijkstra( T src, T dest) {
    if (nodes.find(src) == nodes.end()) {
        throw invalid_argument("Source node does not exist");
    }
    if (nodes.find(dest) == nodes.end()) {
        throw invalid_argument("Destination node does not exist");
    }
    if (src == dest) {
        throw invalid_argument("Source and destination are the same");
    }
    vector<double> distances(V+deletedVertices, INT_MAX);
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

            if (hasEdge(curNode, neighbor)) {
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
    }
    if (distances[dest] == INT_MAX) {
        return -1;
    }
    return distances[dest];
}
// ________________________end dijkstra_________________________________________
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
    // check all the weights to make sure no negative weights
    for (auto node : nodes) {
        for (auto neighbor : node.second->neighbors) {
            if (neighbor.second < 0) {
                cout << "Negative edge weights are not allowed" << endl;
                return vector<pair<T,T>>();
            }
        }
    }
    vector<double> distances(V+deletedVertices, INT_MAX);
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

            if (hasEdge(curNode, neighbor)) {
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
    }
    if (print) {
        for ( auto p : path) {
            cout << p.first << " -> " << p.second << " = " << distances[p.second] << endl;
        }
    }
    return path;
}
// ________________________end dijkstra_________________________________________
/**
 * @brief finds the shortest path between two nodes, if it exists
 * @param src  : source node
 * @param dest  : destination node
 * @param print  : true if the path should be printed, false otherwise
 * @return  : vector of nodes in the shortest path, empty vector if no path exists
 */
template<class T> 
vector<T> Graph<T>::shortestPath(T src, T dest, bool print) {
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
}
// ____________________________end shortestPath_________________________________
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
    if (weighted) {
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
}
// ____________________________end shortestPaths________________________________
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
}
// ______________________________end getComponents______________________________


//##############################################################################
//***************************PRIVATE METHODS************************************
//##############################################################################

/**
 * @brief calculates the average degree of the graph which is (2*E)/V
 * @tparam T
 * @return
 */
template<class T>
double Graph<T>::calculateAverageDegree() {
    double sum = 0;
    // calculate the average degree for when it is directed or not directed
    for (auto node : nodes) {
        sum += node.second->neighbors.size();
    }
    if (isDirected) {
        return (2 * sum) / V;
    } else {
        return (sum) / V;
    }
}
// ________________________end calculateAverageDegree___________________________
/**
 * @brief calculates the average path length of the graph which is (E)/V
 * @tparam T
 * @return
 */
template<class T>
double Graph<T>::calculateAveragePathLength() {
    double sum = 0;
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < V; j++) {
            auto pl = pathLength(i, j);
            if (pl != -1) {
                sum += pl;
            }
        }
    }
    return sum / (V * (V - 1));

}
// ________________________end calculateAveragePathLength_______________________
/**
 * @brief calculates the clustering coefficient of a node which is the number
 *  of triangles divided by the number of possible triangles or V(V-1)/2
 * @param node the node to calculate the clustering coefficient of
 * @return the clustering coefficient of the node
 */
template<class T>
double Graph<T>::clusteringCoef(std::shared_ptr<GraphNode<T>>  node) {
    // use hasEdge and the number of neighbors to calculate the clustering
    // coefficient
    if (node->neighbors.size() < 2) {
        return 0;
    }
    // print out degree of node
    int degree = node->neighbors.size();
    int possibleTriangles = node->neighbors.size() * (node->neighbors.size() - 1);
    int triangles = 0;
    for (auto neighbor : node->neighbors) {
        for (auto neighbor2 : node->neighbors) {
            if (neighbor.first != neighbor2.first && hasEdge(neighbor.first, neighbor2.first)) {
                triangles++;
            }
        }
    }
    return (double) triangles / possibleTriangles;
}
// ________________________end clusteringCoef___________________________________
/**
 * @brief calculates the average global clustering coefficient of the graph
 * @return average global clustering coefficient
 */
template<class T>
double Graph<T>::calculateGlobalClusteringCoefficient() {
    double sum = 0;
    for ( auto node : nodes) {
        sum += clusteringCoef(node.second);
    }
    return sum / V;
}
// ________________end calculateGlobalClusteringCoefficient_____________________

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
}
// ________________________end bfsUtil__________________________________________
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
}
// ___________________________end dfsUtil_______________________________________
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
}
// __________________________end dfsUtil________________________________________
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
                             bool *visited,
                             bool *recur,
                             T parent) {
    if (isDirected) {
        if (!visited[node]) {
            visited[node] = true;
            recur[node] = true;
            // recur all the vertices adjacent to this vertex
            for ( auto neighbor: nodes[node]->neighbors) {
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
        for ( auto neighbor: nodes[node]->neighbors) {
            if (!visited[neighbor.first]) {
                 bool cycle = hasCycleUtil(neighbor.first,
                                                           visited, recur,
                                                           node);
                if (cycle) {
                    return true;
                }
            } else if (neighbor.first != parent) {
                return true;
            }
        }
        return false;
    }
}
// ________________________end isCyclicUtil_____________________________________
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
}
// ________________________end isConnectedUtil__________________________________
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
}
// ________________________end isConnectedUtil__________________________________

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
}
// ________________________end traverse_________________________________________

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
}
// ________________________end cycleFromVertexUtil______________________________
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
}
// ________________________end findBackEdgesUtil________________________________
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
}
// ________________________end getNumberOfComponents____________________________
template<class T>
int Graph<T>::getDegree(T node) {
    auto it = nodes.find(node);
    if (it != nodes.end()) {
        return nodes[node]->neighbors.size();
    }
    return -1;
}
// ______________________________end degree_____________________________________
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
}
// _______________________________end getId_____________________________________
/**
 * @brief method to determine if the graph is bipartite
 * @param node : node to start the search from
 * @param visited : vector of visited nodes
 * @param color  : color array to store the color of each node
 * @return  true if the graph is bipartite, false otherwise
 */
template<class T>
bool Graph<T>::isBipartiteUtil( T node,
                                bool *visited,
                                int *color) {
    visited[node] = true;
    color[node] = 1;
    // checking if a directed graph is bipartite
    for ( auto neighbor : nodes[node]->neighbors) {
        if (!visited[neighbor.first]) {
            if (!isBipartiteUtil(neighbor.first, visited, color)) {
                return false;
            }
        } else if (color[neighbor.first] == color[node]) {
            return false;
        }
    }
    return true;
}
// __________________________end isBipartiteUtil________________________________
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
// __________________________end encode_________________________________________
template<class T>
double Graph<T>::pathLengthUtil(T node1, T node2) {
    if (node1 == node2) {
        return 0;
    }
    std::unique_ptr<bool[]> visited = std::make_unique<bool[]>(V+ deletedVertices);
    std::unique_ptr<double[]> distance = std::make_unique<double[]>(V+ deletedVertices);
    for (int i = 0; i < V+deletedVertices; i++) {
        distance[i] = INT_MAX;
        visited[i] = false;
    }
    // print out the visited and distance array to check its initial values


    if (isWeighted) {
        distance[node1] = 0;
        queue<T> q;
        q.push(node1);
        while (!q.empty()) {
            auto u = q.front();
            q.pop();
            for ( auto neighbor : nodes[u]->neighbors) {
                // make sure the two vertexes exist in the graph
                if (hasEdge(u, neighbor.first)) {
                    if (!visited[neighbor.first]) {
                        visited[neighbor.first] = true;
                        distance[neighbor.first] = distance[u] + neighbor.second;
                        q.push(neighbor.first);
                    }
                }
            }
        }
    } else {
        distance[node1] = 0;
        queue<T> q;
        q.push(node1);
        while (!q.empty()) {
            auto u = q.front();
            q.pop();
            for ( auto neighbor : nodes[u]->neighbors) {
                if (!visited[neighbor.first] && distance[neighbor.first] > distance[u] + 1) {
                    distance[neighbor.first] = distance[u] + 1;
                    q.push(neighbor.first);
                }
            }
            visited[u] = true;
        }
    }
    auto length = distance[node2];
    if (length == INT_MAX || length < 0) {
        length = -1;
    }
    return length;
}
// __________________________end pathLengthUtil_________________________________
template<class T>
void Graph<T>::iteratorPrinter() {
    // use iterator to print the graph
    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        cout << it->first << ": ";
        for (auto it2 = it->second->neighbors.begin(); it2 !=
                  it->second->neighbors.end(); it2++) {
            cout << it2->first << " ";
        }
        cout << endl;
    }
}
// __________________________end iteratorTest___________________________________
template<class T>
Graph<T>::~Graph() {
    // make sure all the smart pointers are deleted
    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        it->second.reset();
    }
}

template<class T>
vector<T> Graph<T>::topologicalSort() {
    // check that the graph is a directed acyclic graph
    if (!isDAG()) {
        cout << "Graph is not a directed acyclic graph" << endl;
        return vector<T>();
    }
    // use a stack to store the nodes in topological order
    vector<T> topologicalOrder;
    std::unique_ptr<bool[]> visited = std::make_unique<bool[]>(V+ deletedVertices);
    std::unique_ptr<int[]> inDegree = std::make_unique<int[]>(V+ deletedVertices);
    for (int i = 0; i < V+deletedVertices; i++) {
        inDegree[i] = 0;
        visited[i] = false;
    }
    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        for (auto it2 = it->second->neighbors.begin(); it2 !=
                  it->second->neighbors.end(); it2++) {
            inDegree[it2->first]++;
        }
    }
    queue<T> q;
    for (int i = 0; i < V+deletedVertices; i++) {
        if (inDegree[i] == 0) {
            q.push(i);
        }
    }
    while (!q.empty()) {
        auto u = q.front();
        q.pop();
        visited[u] = true;
        topologicalOrder.push_back(u);
        for (auto it = nodes[u]->neighbors.begin(); it !=
                  nodes[u]->neighbors.end(); it++) {
            inDegree[it->first]--;
            if (inDegree[it->first] == 0) {
                q.push(it->first);
            }
        }
    }
    return topologicalOrder;
}

template<class T>
vector<vector<T>> Graph<T>::stronglyConnectedComponents() {
    // use Tarjan's algorithm to find the strongly connected components
    vector<vector<T>> scc;
    int* disc = new int[V+deletedVertices];
    int* low = new int[V+deletedVertices];
    bool* stackMember = new bool[V+deletedVertices];
    for (int i = 0; i < V+deletedVertices; i++) {
        disc[i] = -1;
        low[i] = -1;
        stackMember[i] = false;
    }
    stack<T>* st = new stack<T>();
    for (int i = 0; i < V+deletedVertices; i++) {
        if (disc[i] == -1) {
            sccUtil(i, disc, low,  st, stackMember, scc);
        }
    }
    delete[] disc;
    delete[] low;
    delete[] stackMember;
    delete st;
    return scc;
}

template<class T>
Graph<T> Graph<T>::getTranspose() {
    // make and adjacency list from the graph then transpose it
    auto adj = getAdjMatrix();
    Graph<T> g(V+deletedVertices, this->isDirected);
    for (int i = 0; i < adj.size(); i++) {
        for (int j = 0; j < adj[i].size(); j++) {
            if (adj[i][j] == 1) {
                g.addEdge(j, i);
            }
        }
    }
    return g;
}

template<class T>
void Graph<T>::fillOrder(int v, bool *visited, stack<T> &Stack) {
    // use DFS to fill the stack with the nodes in topological order
    visited[v] = true;
    for (auto it = nodes[v]->neighbors.begin(); it !=
              nodes[v]->neighbors.end(); it++) {
        if (!visited[it->first]) {
            fillOrder(it->first, visited, Stack);
        }
    }
    Stack.push(v);
}

template<class T>
void Graph<T>::dfsUtil(int v, bool *visited, vector<T> &comp) {
    visited[v] = true;
    comp.push_back(v);
    for (auto it = nodes[v]->neighbors.begin(); it !=
              nodes[v]->neighbors.end(); it++) {
        if (!visited[it->first]) {
            dfsUtil(it->first, visited, comp);
        }
    }
}

template<class T>
void Graph<T>::addNode(GraphNode<T> *node) {
    // add a node to the graph
    nodes[node->id] = node;
}

template<class T>
vector<vector<T>> Graph<T>::getAdjMatrix() {
    // return the adjacency matrix
    vector<vector<T>> adj(V+deletedVertices, vector<T>(V+deletedVertices, 0));
    for (auto it = nodes.begin(); it != nodes.end(); it++) {
        for (auto it2 = it->second->neighbors.begin(); it2 !=
                  it->second->neighbors.end(); it2++) {
            if (isWeighted) {
                adj[it->first][it2->first] = it2->second;
            } else {
                adj[it->first][it2->first] = 1;
            }
        }
    }
    return adj;
}

template<class T>
vector<pair<T, T>> Graph<T>::bellmanFord(T src, bool print) {
    // use Bellman-Ford algorithm to find the shortest path from src to dest
    vector<pair<T, T>> shortestPath;
    int v = V+deletedVertices;
    int e = E;
    vector<pair<pair<T,T>,double>> edges = getEdges();
    vector<double> dist(v, numeric_limits<double>::max());
    dist[src] = 0;
    for (int i = 0; i < v-1; i++) {
        for (int j = 0; j < e; j++) {
            int s = edges[j].first.first;
            int t = edges[j].first.second;
            double weight = edges[j].second;
            if (dist[s] != numeric_limits<double>::max() &&
                dist[s] + weight < dist[t]) {
                dist[t] = dist[s] + weight;
            }
        }
    }
    for (int j = 0; j < e; j++) {
        int s = edges[j].first.first;
        int t = edges[j].first.second;
        double weight = edges[j].second;
        if (dist[s] != numeric_limits<double>::max() &&
            dist[s] + weight < dist[t]) {
            if (print) {
                cout << "Graph contains negative weight cycle" << endl;
            }
            return shortestPath;
        }
    }
    for (int i = 0; i < v; i++) {
        if (dist[i] != numeric_limits<double>::max()) {
            shortestPath.push_back(make_pair(i, dist[i]));
        }
    }
    if (print) {
        printf("Vertex -- Distance from %d:\n", src);
        for (int i = 0; i < v; i++) {
            if (dist[i] != numeric_limits<double>::max()) {
                printf("  %d   \t      %f\n", i, dist[i]);
            }
        }
    }
    return shortestPath;
}

template<class T>
vector<pair<T, T>> Graph<T>::getBridges() {
    vector<pair<T, T>> bridges;
    bool *visited = new bool[V+deletedVertices];
    int *disc = new int[V+deletedVertices];
    int *low = new int[V+deletedVertices];
    int *parent = new int[V+deletedVertices];
    for (int i = 0; i < V+deletedVertices; i++) {
        visited[i] = false;
        disc[i] = numeric_limits<int>::max();
        low[i] = numeric_limits<int>::max();
        parent[i] = -1;
    }
    for (int i = 0; i < V+deletedVertices; i++) {
        if (!visited[i]) {
            bridgeUtil(i, visited, disc, low, parent, bridges);
        }
    }
    delete[] visited;
    delete[] disc;
    delete[] low;
    delete[] parent;
    return bridges;
}

template<class T>
void
Graph<T>::sccUtil(int u, int *disc, int *low, stack<T> *st, bool *stackMember,
                  vector<vector<T>> &scc) {
    static int time = 0;
    disc[u] = low[u] = ++time;
    st->push(u);
    stackMember[u] = true;
    for (auto it = nodes[u]->neighbors.begin(); it !=
              nodes[u]->neighbors.end(); it++) {
        if (disc[it->first] == -1) {
            sccUtil(it->first, disc, low, st, stackMember, scc);
            low[u] = min(low[u], low[it->first]);
        } else if (stackMember[it->first]) {
            low[u] = min(low[u], disc[it->first]);
        }
    }
    if (low[u] == disc[u]) {
        vector<T> scc_temp;
        while (1) {
            T x = st->top();
            st->pop();
            stackMember[x] = false;
            scc_temp.push_back(x);
            if (x == u) {
                break;
            }
        }
        scc.push_back(scc_temp);
    }
}

template<class T>
void
Graph<T>::bridgeUtil(int u, bool *visited, int *disc, int *low, int *parent,
                     vector<pair<T, T>> &bridges) {
    static int time = 0;
    visited[u] = true;
    disc[u] = low[u] = ++time;
    for (auto it = nodes[u]->neighbors.begin(); it !=
              nodes[u]->neighbors.end(); it++) {
        if (!visited[it->first]) {
            parent[it->first] = u;
            bridgeUtil(it->first, visited, disc, low, parent, bridges);
            low[u] = min(low[u], low[it->first]);
            if (low[it->first] > disc[u]) {
                bridges.push_back(make_pair(u, it->first));
            }
        } else if (it->first != parent[u]) {
            low[u] = min(low[u], disc[it->first]);
        }
    }
}

template<class T>
vector<T> Graph<T>::getArticulationPoints() {
    vector<T> articulationPoints;
    bool *visited = new bool[V+deletedVertices];
    int *disc = new int[V+deletedVertices];
    int *low = new int[V+deletedVertices];
    int *parent = new int[V+deletedVertices];
    for (int i = 0; i < V+deletedVertices; i++) {
        visited[i] = false;
        disc[i] = numeric_limits<int>::max();
        low[i] = numeric_limits<int>::max();
        parent[i] = -1;
    }
    for (int i = 0; i < V+deletedVertices; i++) {
        if (!visited[i]) {
            articulationPointUtil(i, visited, disc, low, parent,
                                  articulationPoints);
        }
    }
    delete[] visited;
    delete[] disc;
    delete[] low;
    delete[] parent;
    return articulationPoints;
}

template<class T>
void Graph<T>::articulationPointUtil(int i, bool *visited, int *disc, int *low,
                                     int *parent, vector<T> &articulationPoints) {
    static int time = 0;
    visited[i] = true;
    disc[i] = low[i] = ++time;
    int children = 0;
    for (auto it = nodes[i]->neighbors.begin(); it !=
              nodes[i]->neighbors.end(); it++) {
        if (!visited[it->first]) {
            children++;
            parent[it->first] = i;
            articulationPointUtil(it->first, visited, disc, low, parent,
                                  articulationPoints);
            low[i] = min(low[i], low[it->first]);
            if (parent[i] == -1 && children > 1) {
                articulationPoints.push_back(i);
            } else if (parent[i] != -1 && low[it->first] >= disc[i]) {
                articulationPoints.push_back(i);
            }
        } else if (it->first != parent[i]) {
            low[i] = min(low[i], disc[it->first]);
        }
    }
}


#endif //GRAPH_CPP_GRAPH_H
// ____________________________Static functions_________________________________


void dfs_utl(vector<vector<int>> mat, vector<vector<bool>>& visited,
                       vector<vector<int>>& mem, int i, int j, int m, int n) {

    visited[i][j] = true;

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

bool dfs_utl(vector<vector<int>>& graph, vector<int> colors, int i, int c) {
    if(colors[i]!=-1){
        return colors[i]==c;
    }
    colors[i] = c;
    for(int j=0;j<graph[i].size();j++){
        if(!dfs_utl(graph,colors,graph[i][j],1-c)){
            return false;
        }
    }
    return true;
}
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
    return g.cycle();
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
    auto size = n * m;
    Graph<int> graph(grid, size, true);
    int total = graph.dijkstra(0, size - 1) + grid[0][0];
    return total;
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

int longestPathSequence(vector<vector<int> > matrix){
    //return the longest path sequence
    auto m = matrix.size();
    auto n = matrix[0].size();
    vector<vector<bool> > visited(m+1,vector<bool>(n+1,false));
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



