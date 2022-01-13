#ifndef GRAPH_CPP_MAXFLOWGRAPH_H
#define GRAPH_CPP_MAXFLOWGRAPH_H
#include <bits/stdc++.h>
using namespace std;
#define FLOATING_POINT_EPSILON 1e-11
    /**
     *  The FlowEdge class represents a capacitated edge with a
      * flow in a FlowNetwork. Each edge consists of two integers
     *  (naming the two vertices), a real-valued capacity, and a real-valued
     *  flow. The data type provides methods for accessing the two endpoints
     *  of the directed edge and the weight. It also provides methods for
     *  changing the amount of flow on the edge and determining the residual
     *  capacity of the edge.
     *  <p>
     *  For additional documentation, see
     *  <a href="https://algs4.cs.princeton.edu/64maxflow">Section 6.4</a> of
     *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
     *
     *  @author Robert Sedgewick and Kevin Wayne (Java)
     *  @author re-implemented in C++ by Ryan Zurrin
     */
    class FlowEdge {
    private:
        int v{}; // edge source
        int w{}; // edge target
        double capacity{}; // capacity of edge
        double flow{}; // current flow
    public:
        FlowEdge(int v, int w, const double capacity);
        FlowEdge(int v, int w, const double capacity, const double flow);
        FlowEdge(const FlowEdge &other);
        FlowEdge(FlowEdge &&other);
        FlowEdge &operator=(const FlowEdge &other);
        FlowEdge &operator=(FlowEdge &&other);
        int from() const;
        int to() const;
        double getCapacity() const;
        double getFlow() const;
        void addFlow(double flow_);
        int other(int vertex) const;
        double residualCapacityTo(int vertex) const;
        void addResidualFlowTo(int vertex, double delta);
        string toString() const;
        friend ostream &operator<<(ostream &os, const FlowEdge &e);
    };
    /**
     * @brief Initializes an edge from vertex v to vertex w with
     * the given capacity and zero flow.
     * @param v the tail vertex
     * @param w the head vertex
     * @param capacity the capacity of the edge
     * @throws IllegalArgumentException if either v or w is a negative integer
     * @throws IllegalArgumentException if capacity is negative
     */
    FlowEdge::FlowEdge(int v, int w, const double capacity) :
            v(v), w(w), capacity(capacity), flow(0.0) {
        if (v < 0) throw invalid_argument("vertex index must be a nonnegative integer");
        if (w < 0) throw invalid_argument("vertex index must be a nonnegative integer");
        if (capacity < 0) throw invalid_argument("capacity must be nonnegative");
    }
    FlowEdge::FlowEdge(int v, int w, const double capacity, const double flow) :
            v(v), w(w), capacity(capacity), flow(flow) {
        if (v < 0) throw std::invalid_argument("vertex < 0");
        if (w < 0) throw std::invalid_argument("vertex < 0");
        if (!(capacity >= 0.0))throw std::invalid_argument("capacity < 0");
        if (!(flow <= capacity)) throw std::invalid_argument("flow > capacity");
        if (!(flow >= 0.0)) throw std::invalid_argument("flow < 0");
    }
    FlowEdge::FlowEdge(const FlowEdge &other) :
            v(other.v), w(other.w), capacity(other.capacity),
            flow(other.flow) {}
    FlowEdge::FlowEdge(FlowEdge &&other) :
            v(other.v), w(other.w), capacity(other.capacity),
            flow(other.flow) {}
    FlowEdge &FlowEdge::operator=(const FlowEdge &other) {
        if (this != &other) {
            v = other.v;
            w = other.w;
            capacity = other.capacity;
            flow = other.flow;
        }
        return *this;
    }
    FlowEdge &FlowEdge::operator=(FlowEdge &&other) {
        if (this != &other) {
            v = other.v;
            w = other.w;
            capacity = other.capacity;
            flow = other.flow;
        }
        return *this;
    }
    /**
     * @brief Returns the tail vertex of the edge.
     * @return the tail vertex of the edge
     */
    int FlowEdge::from() const {
        return v;
    }
    /**
     * Returns the head vertex of the edge.
     * @return the head vertex of the edge
     */
    int FlowEdge::to() const {
        return w;
    }
    /**
     * @brief Returns the capacity of the edge.
     * @return the capacity of the edge
     */
    double FlowEdge::getCapacity() const {
        return capacity;
    }
    /**
     * @brief Returns the flow on the edge.
     * @return the flow on the edge
     */
    double FlowEdge::getFlow() const {
        return flow;
    }
    void FlowEdge::addFlow(double flow_) {
        this->flow += flow_;
    }
    /**
     * @brief Returns the endpoint of the edge that is different from the given
     * vertex (unless the edge represents a self-loop in which case it returns
     * the same vertex).
     * @param vertex one endpoint of the edge
     * @return the endpoint of the edge that is different from the given vertex
     * @throws IllegalArgumentException if vertex is not one of the endpoints
     *   of the edge
     */
    int FlowEdge::other(int vertex) const {
        if (vertex == v) return w;
        else if (vertex == w) return v;
        else throw runtime_error("invalid vertex");
    }
    /**
     * @brief Returns the residual capacity of the edge in the direction
     *  to the given vertex.
     * @param vertex one endpoint of the edge
     * @return the residual capacity of the edge in the direction to the given vertex
     *   If vertex is the tail vertex, the residual capacity equals
     *   capacity() - flow(); if vertex is the head vertex, the
     *   residual capacity equals flow().
     * @throws IllegalArgumentException if {@code vertex} is not one of the endpoints of the edge
     */
    double FlowEdge::residualCapacityTo(int vertex) const {
        if (vertex == v) return flow;
        else if (vertex == w) return capacity - flow;
        else throw runtime_error("invalid vertex");
    }
    /**
     * @brief Increases the flow on the edge in the direction to the given vertex.
     *   If vertex is the tail vertex, this increases the flow on the edge by delta
     *   if vertex is the head vertex, this decreases the flow on the edge by delta.
     * @param vertex one endpoint of the edge
     * @param delta amount by which to increase flow
     * @throws IllegalArgumentException if vertex is not one of the endpoints
     *   of the edge
     * @throws IllegalArgumentException if  delta makes the flow on
     *   on the edge either negative or larger than its capacity
     * @throws IllegalArgumentException if  delta is  NaN
     */
    void FlowEdge::addResidualFlowTo(int vertex, double delta) {
        if (!(delta >= 0.0)) throw std::invalid_argument("delta < 0");
        if (vertex == v) flow -= delta;
        else if (vertex == w) flow += delta;
        else throw runtime_error("invalid vertex");
        // round flow to 0 or capacity if within floating point error
        if (abs(flow) <= FLOATING_POINT_EPSILON) flow = 0.0;
        if (abs(flow - capacity) <= FLOATING_POINT_EPSILON) flow = capacity;
        if (!(flow >= 0.0)) throw std::invalid_argument("flow < 0");
        if (!(flow <= capacity)) throw std::invalid_argument("flow > capacity");
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



    /**
     *  The {@code FlowNetwork} class represents a capacitated network
     *  with vertices named 0 through <em>V</em> - 1, where each directed
     *  edge is of type {@link FlowEdge} and has a real-valued capacity
     *  and flow.
     *  It supports the following two primary operations: add an edge to the network,
     *  iterate over all of the edges incident to or from a vertex. It also provides
     *  methods for returning the number of vertices <em>V</em> and the number
     *  of edges <em>E</em>. Parallel edges and self-loops are permitted.
     *  <p>
     *  This implementation uses an adjacency-lists representation, which
     *  is a vertex-indexed array of {@link Bag} objects.
     *  All operations take constant time (in the worst case) except
     *  iterating over the edges incident to a given vertex, which takes
     *  time proportional to the number of such edges.
     *  <p>
     *  For additional documentation,
     *  see <a href="https://algs4.cs.princeton.edu/64maxflow">Section 6.4</a> of
     *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
     *
     *  @author Robert Sedgewick and Kevin Wayne (Java)
     *  @author re-implemented in C++ by Ryan Zurrin
     */
    class FlowNetwork {
    private:
        int V{}; // number of vertices
        int E{}; // number of edges
        std::vector<std::vector<FlowEdge>> adj; // adjacency list
    public:
        FlowNetwork(int V);
        FlowNetwork(int V, int E);
        FlowNetwork(const FlowNetwork &network);
        FlowNetwork(FlowNetwork &&network);
        FlowNetwork &operator=(const FlowNetwork &network);
        FlowNetwork &operator=(FlowNetwork &&network);
        // constructor that takes in data from an input stream
        FlowNetwork(int V, istream &in);
        void addEdge(const FlowEdge &e);
        int getV() const;
        int getE() const;
        void validateVertex(int v) const;
        std::vector<FlowEdge> getAdj(int v) const;
        std::vector<FlowEdge> getEdges() const;
        string toString() const;
        friend ostream &operator<<(ostream &os, const FlowNetwork &g);
    };
    /**
     * @brief Initializes an empty flow network with V vertices and 0 edges.
     * @param V the number of vertices
     * @throws IllegalArgumentException if V < 0
     */
    FlowNetwork::FlowNetwork(int V) : V(V) {
        if (V < 0)
            throw runtime_error("number of vertices must be nonnegative");
        adj.resize(V);
        this->E = 0;
        for (int v = 0; v < V; v++) {
            adj[v] = std::vector<FlowEdge>();
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
        if (V < 0)
            throw runtime_error("number of vertices must be nonnegative");
        if (E < 0) throw runtime_error("number of edges must be nonnegative");
        this->E = E;
        adj.resize(V);
        for (int v = 0; v < V; v++) {
            adj[v] = std::vector<FlowEdge>();
        }
        // add random edges using the seeded random number generator
        srand(time(NULL));
        for (int i = 0; i < E; i++) {
            int v = rand() % V;
            int w = rand() % V;
            int capacity = rand() % 100;
            addEdge(FlowEdge(v, w, capacity));
        }
    }
    FlowNetwork::FlowNetwork(const FlowNetwork &network) :
    V(network.V), E(network.E), adj(network.adj) {}
    FlowNetwork::FlowNetwork(FlowNetwork &&network) :
    V(network.V), E(network.E), adj(std::move(network.adj)) {}
    FlowNetwork &FlowNetwork::operator=(const FlowNetwork &network) {
        if (this != &network) {
            V = network.V;
            E = network.E;
            adj = network.adj;
        }
        return *this;
    }
    FlowNetwork &FlowNetwork::operator=(FlowNetwork &&network) {
        if (this != &network) {
            V = network.V;
            E = network.E;
            adj = std::move(network.adj);
        }
        return *this;
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
        if (V < 0)
            throw runtime_error("number of vertices must be nonnegative");
        adj.resize(V);
        for (int v = 0; v < V; v++) {
            adj[v] = std::vector<FlowEdge>();
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
    std::vector<FlowEdge> FlowNetwork::getAdj(int v) const {
        validateVertex(v);
        return adj[v];
    }

    /**
     * @brief return list of all edges - excludes self loops
     */
    std::vector<FlowEdge> FlowNetwork::getEdges() const {
        std::vector<FlowEdge> list;
        for (int v = 0; v < V; v++) {
            int selfLoops = 0;
            for (FlowEdge e: adj[v]) {
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
            for (FlowEdge e: adj[v]) {
                ss << e << " ";
            }
            ss << std::endl;
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



    /**
     *  @brief The FordFulkerson class represents a data type for computing a
     *  maximum st-flow and minimum st-cut in a flow network.
     *  <p>
     *  This implementation uses the <em>Ford-Fulkerson</em> algorithm with
     *  the <em>shortest augmenting path</em> heuristic.
     *  The constructor takes <em>O</em>(<em>E V</em> (<em>E</em> + <em>V</em>))
     *  time, where <em>V</em> is the number of vertices and <em>E</em> is
     *  the number of edges. In practice, the algorithm will run much faster.
     *  The inCut() and value() methods take &Theta(1) time.
     *  It uses &Theta(V) extra space (not including the network).
     *  <p>
     *  This correctly computes the maxflow and mincut if all arithmetic
     *  performed is without floating-point rounding error or arithmetic
     *  overflow. This is guaranteed to be the case if all edge capacities
     *  and initial flow values are integers and the value of the maxflow
     *  does not exceeds 2<sup>52</sup>.
     *  <p>
     *  For additional documentation, see
     *  <a href="https://algs4.cs.princeton.edu/64maxflow">Section 6.4</a> of
     *  <i>Algorithms, 4th Edition</i> by Robert Sedgewick and Kevin Wayne.
     *
     *  @author Robert Sedgewick
     *  @author Kevin Wayne
     *  @author re-implemented in C++ by Ryan Zurrin
     */
    class FordFulkerson {
    private:
        const int V; // number of vertices
        const FlowNetwork &G;
        std::vector<double> excess;
        std::vector<int> height;
        std::vector<int> active;
        std::vector<int> current;
        std::vector<FlowEdge> edgeTo;
        std::vector<bool> marked;
        double value;
    public:
        FordFulkerson(const FlowNetwork &G);
        FordFulkerson(const FlowNetwork &G, int s, int t);
        bool hasAugmentingPath(int s, int t);
        bool inCut(int v) const;
        void validateVertex(int v) const;
        double maxFlow(int s, int t);
        double getValue() const;
        std::vector<FlowEdge> getMaxFlow();
        void bfs(int s);
        void relabel(int v);
        void discharge(int v);
        void push(FlowEdge &e);
        void gap(int s);
        double getExcess(int v) const;
        bool isFeasible(std::vector<FlowEdge> &maxFlow);
        bool isFeasible(int s, int t);
        bool check(int s, int t);
        void initPreflow(int s);
        void printMaxFlow() const;
        void printGraph() const;
        string toString() const;
        friend ostream &operator<<(ostream &os, const FordFulkerson &ff);

        void printMinCut();
    };

/**
 * @brief Constructor
 * @param G the graph
 */
    FordFulkerson::FordFulkerson(const FlowNetwork &G) :
            V(G.getV()), G(G){
        excess.resize(V);
        height.resize(V);
        active.resize(V);
        current.resize(V);
        marked.resize(V);
        for (int v = 0; v < V; v++) {
            height[v] = 0;
            active[v] = v;
            current[v] = 0;
            marked[v] = false;
        }
        value = 0;
        initPreflow(0);
    }

/**
 * @brief Constructor
 * @param G the flow network
 * @param s source vertex
 * @param t sink vertex
 */
    FordFulkerson::FordFulkerson(const FlowNetwork &G, int s, int t) :
            V(G.getV()), G(G){
        excess.resize(V);
        height.resize(V);
        active.resize(V);
        current.resize(V);
        initPreflow(s);
        if (s == t) {
            cout << "s == t" << endl;
            return;
        }
        if (!isFeasible(s, t)) {
            cout << "Initial flow is infeasible" << endl;
           // return;
        }
        value = getExcess(t);
        while (hasAugmentingPath(s, t)) {
            double bottleneck = MAXFLOAT;
            for (int v = t; v != s; v = edgeTo[v].other(v)) {
                bottleneck = min(bottleneck, edgeTo[v].residualCapacityTo(v));
            }
            // augment flow
            for (int v = t; v != s; v = edgeTo[v].other(v)) {
                edgeTo[v].addResidualFlowTo(v, bottleneck);
            }
            value += bottleneck;
        }
        assert (check(s, t));
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
        edgeTo = std::vector<FlowEdge>(V, FlowEdge(0, 0, 0));
        marked = std::vector<bool>(V, false);
        queue<int> q;
        q.push(s);
        marked[s] = true;
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            for (FlowEdge e: G.getAdj(v)) {
                int w = e.other(v);
                // residual capacity from v to w
                if (e.residualCapacityTo(w) > 0 && !marked[w]) {
                    edgeTo[w] = e;
                    marked[w] = true;
                    q.push(w);
                }
            }
        }
        return marked[t];
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
    std::vector<FlowEdge> FordFulkerson::getMaxFlow() {
        std::vector<FlowEdge> list;
        for (int v = 0; v < V; v++) {
            for (FlowEdge e: G.getAdj(v)) {
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
            for (FlowEdge e: G.getAdj(v)) {
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
        for (FlowEdge e: G.getAdj(v)) {
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
        for (FlowEdge e: G.getAdj(v)) {
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
            } else {
                w = e.other(w);
                if (height[w] == height[v] + 1) {
                    current[v] = w;
                } else {
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
        for (FlowEdge e: G.getAdj(v)) {
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
    bool FordFulkerson::isFeasible(std::vector<FlowEdge> &maxFlow) {
        for (FlowEdge e: maxFlow) {
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
        for (FlowEdge e: G.getAdj(s)) {
            if (e.getFlow() < -FLOATING_POINT_EPSILON ||
                e.getFlow() > e.getCapacity() + FLOATING_POINT_EPSILON)
                std::cout << "Edge does not satisfy capacity constraints: " << e
                     << std::endl;
            return false;
        }
        // check that net flow into a vertex equals zero, except at source and sink
        if (abs(value + getExcess(s)) > FLOATING_POINT_EPSILON) {
            std::cout << "Excess at source = " << getExcess(s) << std::endl;
            std::cout << "Max flow         = " << value << std::endl;
            return false;
        }
        if (abs(value - getExcess(t)) > FLOATING_POINT_EPSILON) {
            std::cout << "Excess at sink   = " << getExcess(t) << std::endl;
            std::cout << "Max flow         = " << value << std::endl;
            return false;
        }
        for (int v = 0; v < V; v++) {
            // check that net flow out of vertex equals zero
            if (v == s || v == t) continue;
            else if (abs(getExcess(v)) > FLOATING_POINT_EPSILON) {
                std::cout << "Net flow out of " << v << " doesn't equal zero"
                     << std::endl;
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
            std::cout << "Flow is infeasible" << std::endl;
            return false;
        }
        // check that s is on the source side of min cut and that t is not on source side
        if (!inCut(s)) {
            std::cout << "source " << s << " is not on source side of min cut"
                 << std::endl;
            return false;
        }
        if (inCut(t)) {
            std::cout << "sink " << t << " is on source side of min cut" << std::endl;
            return false;
        }
        // check that value of min cut = value of max flow
        double mincutValue = 0.0;
        for (int v = 0; v < V; v++) {
            for (FlowEdge e: G.getAdj(v)) {
                if (e.from() == v) {
                    if (inCut(e.from()) && inCut(e.to()))
                        mincutValue += e.getCapacity();
                } else {
                    if (inCut(e.from()) && inCut(e.to()))
                        mincutValue -= e.getCapacity();
                }
            }
        }
        if (abs(mincutValue - value) > FLOATING_POINT_EPSILON) {
            std::cout << "Max flow value = " << value << " but min cut value = "
                 << mincutValue << std::endl;
            return false;
        }
        return true;
    }


    void FordFulkerson::initPreflow(int s) {
        for (int v = 0; v < V; v++) {
            if (v == s) {
                current[v] = G.getAdj(v).size();
            } else {
                current[v] = 0;
            }
        }
    }

    string FordFulkerson::toString() const {
        stringstream ss;
        // print all the network data to the stringstream
        ss << "FordFulkerson" << std::endl;
        ss << "V = " << V << std::endl;
        ss << "E = " << G.getE() << std::endl;
        ss << "value = " << getValue() << std::endl;
        ss << "edges = " << std::endl;
        for (int v = 0; v < V; v++) {
            for (FlowEdge e: G.getAdj(v)) {
                if (e.getFlow() > 0) {
                    ss << e << std::endl;
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
        std::cout << "Max flow = " << getValue() << std::endl;
        std::cout << "Edges in max flow = " << std::endl;
        for (int v = 0; v < V; v++) {
            for (FlowEdge e: G.getAdj(v)) {
                if (v == e.from() && e.getFlow() > 0) {
                    std::cout << " " << e;
                }
            }
        }
    }

    void FordFulkerson::printGraph() const {
        std::cout << G << std::endl;
    }

    void FordFulkerson::printMinCut() {
        std::cout << "Min cut = " << std::endl;
        for (int v = 0; v < V; v++) {
            if (inCut(v)) {
                std::cout << v << " ";
            }
        }
        std::cout << std::endl;
    }

#endif //GRAPH_CPP_MAXFLOWGRAPH_H
