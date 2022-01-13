//
// Created by Ryan.Zurrin001 on 1/12/2022.
//

#ifndef GRAPH_CPP_MAXFLOWGRAPH_H
#define GRAPH_CPP_MAXFLOWGRAPH_H
#include <bits/stdc++.h>
using namespace std;

class FlowEdge {
private:
    const int v{}; // edge source
    const int w{}; // edge target
    const double capacity{}; // capacity of edge
    double flow{}; // current flow
public:
    FlowEdge(int v, int w, const double capacity);
    int from() const;
    int to() const;
    double getCapacity() const;
    double getFlow() const;
    void addFlow(double flow_);
    int other(int vertex) const;
    double residualCapacityTo(int vertex) const;
    void addResidualFlowTo(int vertex, double delta);
    string toString() const;
    friend ostream& operator<<(ostream& os, const FlowEdge& e);
};
FlowEdge::FlowEdge(int v, int w, const double capacity) :
v(v), w(w), capacity(capacity), flow(0.0) {}
int FlowEdge::from() const {
    return v;
}
int FlowEdge::to() const {
    return w;
}
double FlowEdge::getCapacity() const {
    return capacity;
}
double FlowEdge::getFlow() const {
    return flow;
}
void FlowEdge::addFlow(double flow_) {
    this->flow += flow_;
}
int FlowEdge::other(int vertex) const {
    if (vertex == v) return w;
    else if (vertex == w) return v;
    else throw runtime_error("invalid vertex");
}
double FlowEdge::residualCapacityTo(int vertex) const {
    if (vertex == v) return flow;
    else if (vertex == w) return capacity - flow;
    else throw runtime_error("invalid vertex");
}
void FlowEdge::addResidualFlowTo(int vertex, double delta) {
    if (vertex == v) flow -= delta;
    else if (vertex == w) flow += delta;
    else throw runtime_error("invalid vertex");
}
string FlowEdge::toString() const {
    stringstream ss;
    ss << "(" << v << "->" << w << ") " << flow << "/" << capacity;
    return ss.str();
}
ostream &operator<<(ostream &os, const FlowEdge &e) {
    os << e.toString();
    return os;
}


class FlowNetwork {
private:
    const int V{}; // number of vertices
    int E{}; // number of edges
    vector<vector<FlowEdge>> adj; // adjacency list
public:
    FlowNetwork(int V);
    FlowNetwork(int V, int E);
    // constructor that takes in data from an input stream
    FlowNetwork(int V, istream& in);
    void addEdge(const FlowEdge& e);
    int getV() const;
    int getE() const;
    void validateVertex(int v) const;
    vector<FlowEdge> getAdj(int v) const;
    vector<FlowEdge> getEdges() const;
    string toString() const;
    friend ostream& operator<<(ostream& os, const FlowNetwork& g);
};
/**
 * @brief Initializes an empty flow network with V vertices and 0 edges.
 * @param V the number of vertices
 * @throws IllegalArgumentException if V < 0
 */
FlowNetwork::FlowNetwork(int V) : V(V) {
    if (V < 0) throw runtime_error("number of vertices must be nonnegative");
    adj.resize(V);
    this->E = 0;
    for(int v = 0; v < V; v++) {
        adj[v] = vector<FlowEdge>();
    }
}
/**
 * @brief Initializes a random flow network with V vertices and E edges.
 * The capacities are integers between 0 and 99 and the flow values are zero.
 * @param V the number of vertices
 * @param E the number of edges
 * @throws IllegalArgumentException if {@code V < 0}
 * @throws IllegalArgumentException if {@code E < 0}
 */
FlowNetwork::FlowNetwork(int V, int E) : V(V) {
    if (V < 0) throw runtime_error("number of vertices must be nonnegative");
    if (E < 0) throw runtime_error("number of edges must be nonnegative");
    this->E = E;
    adj.resize(V);
    for(int v = 0; v < V; v++) {
        adj[v] = vector<FlowEdge>();
    }
    // use the random_device class to generate a random seed for the random number generator
    default_random_engine gen;
    uniform_int_distribution<int> capacity(0, 99);
    for (int i = 0; i < E; i++) {
        int v = rand() % V;
        int w = rand() % V;
        double c = capacity(gen);
        addEdge(FlowEdge(v, w, c));
    }
}
/**
 * @brief Initializes a flow network of size V from an input stream.
 * The format is the number of edges E, followed by the E pairs of vertices
 * and edge capacities, with each entry separated by whitespace.
 * @param in the input stream
 * @throws IllegalArgumentException if the endpoints of any edge are not in prescribed range
 * @throws IllegalArgumentException if the number of vertices or edges is negative
 */
FlowNetwork::FlowNetwork(int V, istream &in) : V(V), E(0) {
    if (V < 0) throw runtime_error("number of vertices must be nonnegative");
    adj.resize(V);
    for(int v = 0; v < V; v++) {
        adj[v] = vector<FlowEdge>();
    }
    in >> E;
    if (E < 0) throw runtime_error("number of edges must be nonnegative");
    for (int i = 0; i < E; i++) {
        int v, w;
        double capacity;
        in >> v >> w >> capacity;
        validateVertex(v);
        validateVertex(w);
        addEdge(FlowEdge(v, w, capacity));
    }
}
/**
 * @brief Adds the edge e to the network.
 * @param e the edge
 * @throws IllegalArgumentException unless endpoints of edge are between
 *         0 and V-1
 */
void FlowNetwork::addEdge(const FlowEdge &e) {
    int v = e.from();
    int w = e.to();
    validateVertex(v);
    validateVertex(w);
    adj[v].push_back(e);
    adj[w].push_back(e);
    E++;
}
/**
 * Returns the number of vertices in the edge-weighted graph.
 * @return the number of vertices in the edge-weighted graph
 */
int FlowNetwork::getV() const {
    return V;
}
/**
 * Returns the number of edges in the edge-weighted graph.
 * @return the number of edges in the edge-weighted graph
 */
int FlowNetwork::getE() const {
    return E;
}
/**
 * @brief Validates that vertex is between 0 and V-1.
 * @throws IllegalArgumentException unless  0 <= v < V
 * @param v the vertex
 */
void FlowNetwork::validateVertex(int v) const {
    if (v < 0 || v >= V)
        throw runtime_error("vertex " + to_string(v) +
        " is not between 0 and " + to_string(V - 1));
}
/**
 * @brief Returns the edges incident on vertex v (includes both edges pointing to
 * and from  v).
 * @param v the vertex
 * @return the edges incident on vertex v as an Iterable
 * @throws IllegalArgumentException unless 0 <= v < V
 */
vector<FlowEdge> FlowNetwork::getAdj(int v) const {
    validateVertex(v);
    return adj[v];
}
/**
 * @brief return list of all edges - excludes self loops
 */
vector<FlowEdge> FlowNetwork::getEdges() const {
    vector<FlowEdge> list;
    for (int v = 0; v < V; v++) {
        int selfLoops = 0;
        for (FlowEdge e : adj[v]) {
            if (e.to() != v)
                list.push_back(e);
            else if (e.from() != v)
                list.push_back(e);
            else
                selfLoops++;
        }
    }
    return list;
}
/**
 * @brief build a string representation of the graph
 * @return  a string representation of the graph
 */
string FlowNetwork::toString() const {
    stringstream ss;
    for (int v = 0; v < V; v++) {
        ss << v << ": ";
        for (FlowEdge e : adj[v]) {
            ss << e << " ";
        }
        ss << endl;
    }
    return ss.str();
}
/**
 * @brief overload << operator to print the graph to an output stream
 * @param os  the output stream
 * @param g   the graph
 * @return  the output stream
 */
ostream &operator<<(ostream &os, const FlowNetwork &g) {
    os << g.toString();
    return os;
}




class FordFulkerson {
private:
    const double FLOATING_POINT_EPSILON; // = 1E-11;
    const int V; // number of vertices
    const FlowNetwork& G;
    vector<double> excess;
    vector<int> height;
    vector<int> active;
    vector<int> current;
    vector<int> edgeTo;
    double value;
public:
    FordFulkerson(const FlowNetwork& G);
    FordFulkerson(const FlowNetwork& G, int s, int t);
    bool hasAugmentingPath(int s, int t);
    bool inCut(int v) const;
    void validateVertex(int v) const;
    double maxFlow(int s, int t);
    double getValue() const;
    vector<FlowEdge> getMaxFlow();
    void bfs(int s);
    void relabel(int v);
    void discharge(int v);
    void push(FlowEdge &e);
    void gap(int s);
    double getExcess(int v) const;
    bool isFeasible(vector<FlowEdge> &maxFlow);
    bool isFeasible(int s, int t);
    bool check(int s, int t);
    void initPreflow(int s);
    void printMaxFlow() const;
    void printGraph() const;
    string toString() const;
    friend ostream& operator<<(ostream& os, const FordFulkerson& ff);
};
/**
 * @brief Constructor
 * @param G the graph
 */
FordFulkerson::FordFulkerson(const FlowNetwork &G) :
V(G.getV()), G(G) , FLOATING_POINT_EPSILON(1E-11) {
    excess.resize(V);
    height.resize(V);
    active.resize(V);
    current.resize(V);
    edgeTo.resize(V);
    initPreflow(0);
}
/**
 * @brief Constructor
 * @param G the flow network
 * @param s source vertex
 * @param t sink vertex
 */
FordFulkerson::FordFulkerson(const FlowNetwork &G, int s, int t) :
V(G.getV()), G(G), FLOATING_POINT_EPSILON(1E-11) {
    excess.resize(V);
    height.resize(V);
    active.resize(V);
    current.resize(V);
    edgeTo.resize(V);
    initPreflow(s);
    bfs(s);
    while (excess[t] > 0) {
        if (active[t] < 0) {
            gap(s);
        }
        relabel(t);
        bfs(s);
    }
}
/**
 * @brief Returns true if there is an augmenting path from source s to sink t in
 *      the residual network.
 * @param s  the source vertex
 * @param t  the sink vertex
 * @return  true if there is an augmenting path, false otherwise
 */
bool FordFulkerson::hasAugmentingPath(int s, int t) {
    validateVertex(s);
    validateVertex(t);
    if (s == t)
        return true;
    bfs(s);
    return active[t] != -1;
}
/**
 * @brief Returns true if vertex v is on the s side of the min-cut
 * @param v  the vertex
 * @return  true if vertex v is on the s side of the min-cut
 */
bool FordFulkerson::inCut(int v) const {
    validateVertex(v);
    return height[v] < V;
}
/**
 * @brief Validates that vertex is between 0 and V-1.
 * @throws IllegalArgumentException unless  0 <= v < V
 * @param v the vertex
 */
void FordFulkerson::validateVertex(int v) const {
    if (v < 0 || v >= V)
        throw runtime_error("vertex " + to_string(v) +
        " is not between 0 and " + to_string(V - 1));
}
/**
 * @brief finds the maximum flow from s to t
 * @param s  the source vertex
 * @param t  the sink vertex
 * @return  the maximum flow
 */
double FordFulkerson::maxFlow(int s, int t) {
    validateVertex(s);
    validateVertex(t);
    if (s == t)
        return 0.0;
    initPreflow(s);
    bfs(s);
    while (excess[t] > 0) {
        if (active[t] < 0) {
            gap(s);
        }
        relabel(t);
        bfs(s);
    }
    return value;
}
/**
 * @brief gets the value of the max flow
 * @return  the value of the max flow
 */
double FordFulkerson::getValue() const {
    return value;
}
/**
 * @brief return list of edges in max flow
 * @return  list of edges in max flow
 */
vector<FlowEdge> FordFulkerson::getMaxFlow() {
    vector<FlowEdge> list;
    for (int v = 0; v < V; v++) {
        for (FlowEdge e : G.getAdj(v)) {
            if (e.residualCapacityTo(v) > 0)
                list.push_back(e);
        }
    }
    return list;
}
/**
 * @brief breadth-first search from s
 * @param s  the source vertex
 */
void FordFulkerson::bfs(int s) {
    for (int v = 0; v < V; v++) {
        height[v] = -1;
        active[v] = -1;
    }
    height[s] = 0;
    active[s] = 0;
    queue<int> q;
    q.push(s);
    while (!q.empty()) {
        int v = q.front();
        q.pop();
        for (FlowEdge e : G.getAdj(v)) {
            int w = e.to();
            if (e.residualCapacityTo(w) > 0 && height[w] < 0) {
                height[w] = height[v] + 1;
                active[w] = v;
                q.push(w);
            }
        }
    }
}
/**
 * @brief relabel the vertex v
 * @param v the vertex
 */
void FordFulkerson::relabel(int v) {
    int minHeight = V;
    for (FlowEdge e : G.getAdj(v)) {
        int w = e.to();
        if (e.residualCapacityTo(w) > 0) {
            minHeight = min(minHeight, height[w]);
        }
    }
    height[v] = minHeight + 1;
}
/**
 * @brief Discharge a vertex
 * @param v  Vertex to discharge
 */
void FordFulkerson::discharge(int v) {
    for (FlowEdge e : G.getAdj(v)) {
        int w = e.to();
        if (e.residualCapacityTo(w) > 0) {
            if (height[w] > height[v]) {
                push(e);
                if (excess[v] == 0)
                    return;
            }
        }
    }
    gap(v);
}
/**
 * @brief pushes flow along edge e
 * @param e edge to push flow along
 */
void FordFulkerson::push(FlowEdge &e) {
    int v = e.from();
    int w = e.to();
    double amt = min(excess[v], e.residualCapacityTo(w));
    e.addResidualFlowTo(w, amt);
    excess[v] -= amt;
    excess[w] += amt;
}
/**
 * @brief gap relabels all vertices in the gap
 * @param s  - source
 */
void FordFulkerson::gap(int s) {
    while (true) {
        int v = active[s];
        if (v < 0)
            return;
        int w = current[v];
        FlowEdge &e = reinterpret_cast<FlowEdge &>(edgeTo[v]);
        push(e);
        if (excess[v] == 0) {
            active[s] = -1;
            return;
        }
        else {
            w = e.other(w);
            if (height[w] == height[v] + 1) {
                current[v] = w;
            }
            else {
                active[s] = -1;
            }
        }
    }
}
/**
 * @brief return the excess flow at vertex v
 * @param v  vertex
 * @return  excess flow at vertex v
 */
double FordFulkerson::getExcess(int v) const {
    double excess_ = 0.0;
    for (FlowEdge e : G.getAdj(v)) {
        if (e.from() == v)
            excess_ -= e.getFlow();
        else
            excess_ += e.getFlow();
    }
    return excess_;
}
/**
 * @brief check if the the flow is feasible given the current flow values
 * @param maxFlow  the maximum flow value
 * @return true if the flow is feasible
 */
bool FordFulkerson::isFeasible(vector<FlowEdge> &maxFlow) {
    for (FlowEdge e : maxFlow) {
        if (e.getFlow() < -FLOATING_POINT_EPSILON ||
            e.getFlow() > e.getCapacity() + FLOATING_POINT_EPSILON)
            return false;
    }
    return true;
}
/**
 * @brief check if the given graph is a flow network
 * @param s source vertex
 * @param t sink vertex
 * @return true if the graph is a flow network, false otherwise
 */
bool FordFulkerson::isFeasible(int s, int t) {
    // check that capacity constraints are satisfied
    for (FlowEdge e : G.getAdj(s)) {
        if (e.getFlow() < -FLOATING_POINT_EPSILON ||
            e.getFlow() > e.getCapacity() + FLOATING_POINT_EPSILON)
            cout << "Edge does not satisfy capacity constraints: " << e << endl;
            return false;
    }
    // check that net flow into a vertex equals zero, except at source and sink
    if (abs(value + getExcess(s)) > FLOATING_POINT_EPSILON) {
        cout << "Excess at source = " << getExcess(s) << endl;
        cout << "Max flow         = " << value << endl;
        return false;
    }
    if (abs(value - getExcess(t)) > FLOATING_POINT_EPSILON) {
        cout << "Excess at sink   = " << getExcess(t) << endl;
        cout << "Max flow         = " << value << endl;
        return false;
    }
    for (int v = 0; v < V; v++) {
        // check that net flow out of vertex equals zero
        if (v == s || v == t) continue;
        else if (abs(getExcess(v)) > FLOATING_POINT_EPSILON) {
            cout << "Net flow out of " << v << " doesn't equal zero" << endl;
            return false;
        }
    }
    return true;
}

/**
 * @brief check optimality conditions for max flow
 * @param s  source vertex
 * @param t  sink vertex
 * @return   true if optimality conditions are satisfied, false otherwise
 */
bool FordFulkerson::check(int s, int t) {
    // check that flow is feasible
    if (!isFeasible(s, t)) {
        cout << "Flow is infeasible" << endl;
        return false;
    }
    // check that s is on the source side of min cut and that t is not on source side
    if (!inCut(s)) {
        cout << "source " << s << " is not on source side of min cut" << endl;
        return false;
    }
    if (inCut(t)) {
        cout << "sink " << t << " is on source side of min cut" << endl;
        return false;
    }
    // check that value of min cut = value of max flow
    double mincutValue = 0.0;
    for (int v = 0; v < V; v++) {
        for (FlowEdge e : G.getAdj(v)) {
            if (e.from() == v) {
                if (inCut(e.from()) && inCut(e.to()))
                    mincutValue += e.getCapacity();
            }
            else {
                if (inCut(e.from()) && inCut(e.to()))
                    mincutValue -= e.getCapacity();
            }
        }
    }
    if (abs(mincutValue - value) > FLOATING_POINT_EPSILON) {
        cout << "Max flow value = " << value << " but min cut value = " << mincutValue << endl;
        return false;
    }
    return true;
}


void FordFulkerson::initPreflow(int s) {
    for (int v = 0; v < V; v++) {
        excess[v] = 0;
        height[v] = 0;
        active[v] = -1;
        current[v] = -1;
        edgeTo[v] = -1;
    }
    for (FlowEdge e : G.getAdj(s)) {
        if (e.residualCapacityTo(s) > 0) {
            excess[s] += e.residualCapacityTo(s);
            edgeTo[s] = e.to();
            current[s] = e.to();
        }
    }
    height[s] = V;
    active[s] = -1;
}

string FordFulkerson::toString() const {
    stringstream ss;
    // print all the network data to the stringstream
    ss << "FordFulkerson" << endl;
    ss << "V = " << V << endl;
    ss << "E = " << G.getE() << endl;
    ss << "value = " << getValue() << endl;
    ss << "edges = " << endl;
    for (int v = 0; v < V; v++) {
        for (FlowEdge e : G.getAdj(v)) {
            if (e.getFlow() > 0) {
                ss << e << endl;
            }
        }
    }
    return ss.str();
}
/**
 * @brief  overload the << operator to print the network data to ouput stream
 * @param os  output stream
 * @param ff  FordFulkerson object
 * @return    output stream
 */
ostream &operator<<(ostream &os, const FordFulkerson &ff) {
    os << ff.toString();
    return os;
}
/**
 * @brief print the network data to the console
 */
void FordFulkerson::printMaxFlow() const {
    cout << "Max flow = " << getValue() << endl;
    cout << "Edges in max flow = " << endl;
    for (int v = 0; v < V; v++) {
        for (FlowEdge e : G.getAdj(v)) {
            if (e.getFlow() > 0) {
                cout << e << endl;
            }
        }
    }
}
void FordFulkerson::printGraph() const {
    cout << G << endl;
}


#endif //GRAPH_CPP_MAXFLOWGRAPH_H
