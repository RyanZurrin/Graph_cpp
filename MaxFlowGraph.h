#ifndef GRAPH_CPP_MAXFLOWGRAPH_H
#define GRAPH_CPP_MAXFLOWGRAPH_H
#include <bits/stdc++.h>
#include "Bag.h"
#define DEBUG 0
using namespace std;
#define FLOATING_POINT_EPSILON 1e-10
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
        double flow{};
        // current flow
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
        ~FlowEdge() = default;

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
        if (DEBUG) cout << "residualCapacityTo" << endl;
        if (DEBUG) cout << "v " << v << " w " << w << " vertex " << vertex << endl;
        if (DEBUG) cout << "capacity " << capacity << " flow " << flow << endl;
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
        if (DEBUG) cout << "addResidualFlowTo" << endl;
        if (DEBUG) cout << "delata " << delta << " vertex " << vertex << endl;
        if (vertex == v) flow -= delta;
        if (DEBUG) cout << "flow after -= delta " << flow << endl;
        else if (vertex == w) flow += delta;
        if (DEBUG) cout << "flow after += delta " << flow << endl;
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
        // adjacency list
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
        Bag<FlowEdge *> * getAdj(int v) const;
        Bag<FlowEdge*>* getEdges() const ;
        string toString() const;
        friend ostream &operator<<(ostream &os, const FlowNetwork &g);
        // overload the [] operator to return the Bag of edges
        // incident from a given vertex
        Bag<FlowEdge*>* operator[](int v)  {
            validateVertex(v);
            return this->getAdj(v);
        }
        // overload the () operator to return the Bag of edges
        // incident from a given vertex
        Bag<FlowEdge*>* operator()(int v)  {
            validateVertex(v);
            return this->getAdj(v);
        }

        Bag<FlowEdge*>* adj;
        FlowEdge getEdge(int i, int j) const;

        ~FlowNetwork() {
            if (DEBUG) cout << "destructor called" << endl;
            // iterate over all edges and delete them
            for (int v = 0; v < V; v++) {
                for (Bag<FlowEdge *>::Iterator it = adj[v].begin(); it != adj[v].end(); it++) {
                    delete *it;
                }
            }
            delete[] adj;
        }
    };
    /**
     * @brief Initializes an empty flow network with V vertices and 0 edges.
     * @param V the number of vertices
     * @throws IllegalArgumentException if V < 0
     */
    FlowNetwork::FlowNetwork(int V) : V(V) {
        if (V < 0)
            throw runtime_error("number of vertices must be nonnegative");
        adj = new Bag<FlowEdge*>[V];
        for (int v = 0; v < V; v++) {
            adj[v] = Bag<FlowEdge*>();
        }
    }
    /**
     * @brief Initializes a random flow network with V vertices and E edges.
     * The capacities are integers between 0 and 99 and the flow values are zero.
     * @param V the number of vertices
     * @param E the number of edges
     * @throws IllegalArgumentException if {V < 0}
     * @throws IllegalArgumentException if {E < 0}
     */
    FlowNetwork::FlowNetwork(int V, int E) : V(V) {
        if(DEBUG) cout << "FlowNetwork::FlowNetwork(int V, int E)" << endl;
        if (V < 0)
            throw runtime_error("number of vertices must be nonnegative");
        if (E < 0) throw runtime_error("number of edges must be nonnegative");
        this->E = E;
        if (DEBUG) cout << "V: " << V << " E: " << E << endl;
        adj = new Bag<FlowEdge*>[V];
        for (int v = 0; v < V; v++) {
            if (DEBUG) cout << "v: " << v << endl;
            adj[v] = Bag<FlowEdge*>();
        }
        if (DEBUG) cout << "adj: " << adj << endl;
        // add random edges using the seeded random number generator
        srand(time(NULL));
        for (int i = 0; i < E; i++) {
            if (DEBUG) cout << "i: " << i << endl;
            int v = rand() % V;
            int w = rand() % V;
            int capacity = rand() % 100;
            if (DEBUG) cout << "v: " << v << " w: " << w << " capacity: " << capacity << endl;
            addEdge(FlowEdge(v, w, capacity));
        }
        if (DEBUG) cout << "end of constructor and adj: " << adj << endl;
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
        adj = new Bag<FlowEdge*>[V];
        for (int v = 0; v < V; v++) {
            adj[v] = Bag<FlowEdge*>();
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
            if (DEBUG) cout << "v: " << v << " w: " << w << " capacity: " << capacity << endl;
        }
    }
    /**
     * @brief Adds the edge e to the network.
     * @param e the edge
     * @throws IllegalArgumentException unless endpoints of edge are between
     *         0 and V-1
     */
    void FlowNetwork::addEdge(const FlowEdge &e) {
        if (DEBUG) cout << "FlowNetwork::addEdge(const FlowEdge &e)" << endl;
        int v = e.from();
        int w = e.to();
        validateVertex(v);
        validateVertex(w);
        adj[v].add(new FlowEdge(e));
        adj[w].add(new FlowEdge(e));
        E++;
        if (DEBUG) cout << "end of addEdge and E: " << E << endl;
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
        if (DEBUG) cout << "FlowNetwork::validateVertex(int v)" << endl;
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
    Bag<FlowEdge *> * FlowNetwork::getAdj(int v) const {
        if (DEBUG) cout << "FlowNetwork::getAdj(int v)" << endl;
        validateVertex(v);
        return &adj[v];
    }

    /**
     * @brief return list of all edges - excludes self loops
     */
    Bag<FlowEdge*>* FlowNetwork::getEdges() const {
        if (DEBUG) cout << "FlowNetwork::getEdges()" << endl;
        Bag<FlowEdge*>* edges = new Bag<FlowEdge*>();
        for (int v = 0; v < V; v++) {
            for (FlowEdge* e : *getAdj(v)) {
                if (e->from() != e->to()) {
                    edges->add(e);
                }
            }
        }
        if (DEBUG) cout << "end of getEdges()" << endl;
        return edges;
    }

    /**
     * @brief build a string representation of the graph
     * @return  a string representation of the graph
     */
    string FlowNetwork::toString() const {
        if (DEBUG) cout << "FlowNetwork::toString()" << endl;
        stringstream ss;
        ss << "V = " << V << ", E = " << E << endl;
        for (int v = 0; v < V; v++) {
            ss << "Adj[" << v << "]: ";
            for (FlowEdge* e : *getAdj(v)) {
                ss << e->toString() << "  ";
            }
            ss << endl;
        }
        if (DEBUG) cout << "end of toString()" << endl;
        return ss.str();
    }

    /**
     * @brief overload << operator to print the graph to an output stream
     * @param os  the output stream
     * @param g   the graph
     * @return  the output stream
     */
    ostream &operator<<(ostream &os, const FlowNetwork &g) {
        if (DEBUG) cout << "FlowNetwork::operator<<(ostream &os, const FlowNetwork &g)" << endl;
        os << g.toString();
        return os;
    }

    FlowEdge FlowNetwork::getEdge(int i, int j) const {
        if (DEBUG) cout << "FlowNetwork::getEdge(int i, int j)" << endl;
        for (FlowEdge* e : *getAdj(i)) {
            if (e->to() == j) {
                if (DEBUG) cout << "end of getEdge(int i, int j), returning ";
                if (DEBUG) cout << "e: " << e->toString() << endl;
                return *e;
            }
        }
        if (DEBUG) cout << "end of getEdge(int i, int j), returning ";
        if (DEBUG) cout << "FLOWEDGE(0,0,0)" << endl;
        return FlowEdge(0, 0, 0);
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
        bool* marked;
        vector<FlowEdge*> edgeTo;
        const FlowNetwork& G;
        double value;
    public:
        FordFulkerson(const FlowNetwork &G, int s, int t);
        double getValue() const;
        bool hasAugmentingPath(int s, int t);
        bool inCut(int v) const;
        bool isFeasible(int s, int t);
        void validate(int v) const;
        bool check(int s, int t);
        double getExcess(int v) const;
        void printGraph() const;
        string toString() const;
        friend ostream &operator<<(ostream &os, const FordFulkerson &ff);
        // add iterator to iterate over edges
        class iterator {
        private:
            const FordFulkerson& ff;
            int v;
            int i;
        public:
            iterator(const FordFulkerson& ff, int v, int i) : ff(ff), v(v), i(i) {}
            FlowEdge operator*() const { return ff.G.getEdge(v, i); }
            iterator& operator++() { i++; return *this; }
            bool operator!=(const iterator& other) const { return v != other.v || i != other.i; }
            bool operator==(const iterator& other) const { return v == other.v && i == other.i; }

        }; // end of iterator class
        iterator begin(int v) const {
            return iterator(*this, v, 0);
        }
        iterator end(int v) const {
            return iterator(*this, v, G.getAdj(v)->size());
        }

        ~FordFulkerson() {
            for (FlowEdge* e : *G.getEdges()) {
                delete e;
            }
            delete[] marked;
        }

    };



/**
 * @brief Constructor
 * @param G the flow network
 * @param s source vertex
 * @param t sink vertex
 */
    FordFulkerson::FordFulkerson(const FlowNetwork &G, int s, int t) :
            V(G.getV()), G(G){
        validate(s);
        validate(t);
        if (s == t) throw std::invalid_argument("source equals sink");
        if(!isFeasible(s, t)) throw std::invalid_argument("infeasible");
        value = getExcess(t);
        if (DEBUG) cout << "end of FordFulkerson constructor" << endl;
        if (DEBUG) cout << "value: " << value << endl;
        while (hasAugmentingPath(s, t)) {
            double bottle = std::numeric_limits<double>::max();
            if (DEBUG) cout << "bottle: " << bottle << endl;
            for (int v = t; v != s; v = edgeTo[v]->other(v)) {
                bottle = std::min(bottle, edgeTo[v]->residualCapacityTo(v));
                if (DEBUG) cout << "bottle: " << bottle << endl;
            }
            // augment flow
            for (int v = t; v != s; v = edgeTo[v]->other(v)) {
                if (DEBUG) cout << "adding residualFlowTo: " << edgeTo[v]->residualCapacityTo(v) << endl;
                edgeTo[v]->addResidualFlowTo(v, bottle);
            }
            value += bottle;
            if (DEBUG) cout << "value after += bottle: " << value << endl;
        }
        if (DEBUG) cout << "end of FordFulkerson constructor before assert" << endl;
        assert(check(s, t));
    }

/**
 * @brief Returns true if there is an augmenting path from source s to sink t in
 *      the residual network.
 * @param s  the source vertex
 * @param t  the sink vertex
 * @return  true if there is an augmenting path, false otherwise
 */
    bool FordFulkerson::hasAugmentingPath(int s, int t) {
        if (DEBUG) cout << "hasAugmentingPath(" << s << ", " << t << ")" << endl;
        validate(s);
        validate(t);
        marked = new bool[V];
        edgeTo = vector<FlowEdge*>(V);
        for (int v = 0; v < V; v++) {
            marked[v] = false;
            if (DEBUG) cout << "marked[" << v << "] = " << marked[v] << endl;
        }
        queue<int> q;
        q.push(s);
        marked[s] = true;
        if (DEBUG) cout << "marked[" << s << "] = " << marked[s] << endl;
        while (!q.empty()) {
            int v = q.front();
            q.pop();
            if (DEBUG) cout << "v = " << v << endl;
            for (FlowEdge* e : *G.getAdj(v)) {
                int w = e->other(v);
                if (DEBUG) cout << "w = " << w << endl;
                if (e->residualCapacityTo(w) > 0 && !marked[w]) {
                    if (DEBUG) cout << "marked[" << w << "] = " << marked[w] << endl;
                    if (DEBUG) cout << "edgeTo[" << w << "] = " << edgeTo[w] << endl;
                    edgeTo[w] = e;
                    if (DEBUG) cout << "edgeTo[" << w << "] = " << edgeTo[w] << endl;
                    marked[w] = true;
                    if (DEBUG) cout << "marked[" << w << "] = " << marked[w] << endl;
                    q.push(w);
                    if (DEBUG) cout << "q.push(" << w << ")" << endl;
                }
            }
        }
        if (DEBUG) cout << "marked[" << t << "] = " << marked[t] << endl;
        return marked[t];
    }

/**
 * @brief Returns true if vertex v is on the s side of the min-cut
 * @param v  the vertex
 * @return  true if vertex v is on the s side of the min-cut
 */
    bool FordFulkerson::inCut(int v) const {
        if (DEBUG) cout << "inCut(" << v << ") = " << marked[v] << endl;
        validate(v);
        return marked[v];
    }

/**
 * @brief Validates that vertex is between 0 and V-1.
 * @throws IllegalArgumentException unless  0 <= v < V
 * @param v the vertex
 */
    void FordFulkerson::validate(int v) const {
        if (DEBUG)  cout << "validate(" << v << ")" << endl;
        if (v < 0 || v >= V)
            throw runtime_error("vertex " + to_string(v) +
                                " is not between 0 and " + to_string(V - 1));
    }


/**
 * @brief gets the value of the max flow
 * @return  the value of the max flow
 */
    double FordFulkerson::getValue() const {
        if (DEBUG) cout << "FordFulkerson::getValue() is " << value << endl;
        return value;
    }



/**
 * @brief return the excess flow at vertex v
 * @param v  vertex
 * @return  excess flow at vertex v
 */
    double FordFulkerson::getExcess(int v) const {
        if (DEBUG) cout << "getExcess(" << v << ")" << endl;
        validate(v);
        double excess = 0.0;
        for (FlowEdge* e : *G.getAdj(v)) {
            if (e->from() == v) {
                excess -= e->getFlow();
                if (DEBUG) cout << "excess -= " << e->getFlow() << endl;
            } else {
                excess += e->getFlow();
                if (DEBUG) cout << "excess += " << e->getFlow() << endl;
            }
            if (DEBUG) cout << "excess = " << excess << endl;
        }
        if (DEBUG) cout << "returning excess = " << excess << endl;
        return excess;
    }


/**
 * @brief check if the given graph is a flow network
 * @param s source vertex
 * @param t sink vertex
 * @return true if the graph is a flow network, false otherwise
 */
    bool FordFulkerson::isFeasible(int s, int t) {
        if (DEBUG) {
            cout << "isFeasible" << endl;
            cout << "s: " << s << " t: " << t << endl;
        }
        // check that capacity constraints are satisfied
        for (FlowEdge* e: *G.getAdj(s)) {
            if (e->getFlow() < -FLOATING_POINT_EPSILON ||
                e->getFlow() > e->getCapacity() + FLOATING_POINT_EPSILON) {
                std::cout << "Edge does not satisfy capacity constraints: " << e
                          << std::endl;
                return false;
            }
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
        if (DEBUG) {
            cout << "isFeasible: true" << endl;
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
        if (DEBUG) {
            std::cout << "Checking optimality conditions..." << std::endl;
            std::cout << "Checking flow through " << s << " to " << t << std::endl;
        }
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
            for (FlowEdge* e: *G.getAdj(v)) {
                if (e->from() == v) {
                    if (inCut(e->from()) && inCut(e->to()))
                        mincutValue += e->getCapacity();
                } else {
                    if (inCut(e->from()) && inCut(e->to()))
                        mincutValue -= e->getCapacity();
                }
            }
        }
        if (abs(mincutValue - value) > FLOATING_POINT_EPSILON) {
            std::cout << "Max flow value = " << value << " min cut value = "
                 << mincutValue << std::endl;
            return false;
        }
        if (DEBUG) {
            std::cout << "Checking optimality conditions...OK" << std::endl;
        }
        return true;
    }

    string FordFulkerson::toString() const {
        if (DEBUG) {
            std::cout << "FordFulkerson::toString()" << std::endl;
        }
        stringstream ss;
        // print all the network data to the stringstream
        ss << "Max flow value = " << value << endl;
        ss << "Min cut: " << endl;
        for (int v = 0; v < V; v++) {
            if (inCut(v))
                ss << v << " ";
        }
        ss << endl;
        if (DEBUG) {
            std::cout << "FordFulkerson::toString() done" << std::endl;
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
        if (DEBUG) { os << "FordFulkerson::operator<<" << endl; }
        if (ff.G.getV() == 0) {
            os << "Empty network" << endl;
            return os;
        }
        os << ff.toString();
        if (DEBUG) { os << "FordFulkerson::operator<< end" << endl; }
        return os;
    }


    void FordFulkerson::printGraph() const {
        std::cout << G << std::endl;
    }

#endif //GRAPH_CPP_MAXFLOWGRAPH_H
