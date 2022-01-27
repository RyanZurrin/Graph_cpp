#ifndef GRAPH_CPP_MAXFLOWGRAPH_H
#define GRAPH_CPP_MAXFLOWGRAPH_H
#include <bits/stdc++.h>
#include "Bag.h"
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
    public:
        FlowEdge(int v, int w, double capacity);
        FlowEdge(int v, int w, double capacity, double flow);
        FlowEdge(FlowEdge& e);
        FlowEdge(FlowEdge&& e);
        FlowEdge& operator=(FlowEdge& e);
        FlowEdge& operator=(FlowEdge&& e);
        int from() const;
        int to() const;
        double capacity() const;
        double flow() const;
        int other(int vertex) const;
        double residualCapacityTo(int vertex) const;
        void addResidualFlowTo(int v, double delta);
        string toString() const;
        friend ostream& operator<<(ostream& os, const FlowEdge& e);
    private:
        int _v;
        int _w;
        double _capacity;
        double _flow;
    };

FlowEdge::FlowEdge(int v, int w, double capacity) {
    if (v < 0 || w < 0) {
        throw invalid_argument("vertex index must be nonnegative");
    }
    if (capacity < 0.0) {
        throw invalid_argument("capacity must be nonnegative");
    }
    _v = v;
    _w = w;
    _capacity = capacity;
    _flow = 0.0;
}

FlowEdge::FlowEdge(int v, int w, double capacity, double flow) {
    if (v < 0 || w < 0) {
        throw invalid_argument("vertex index must be nonnegative");
    }
    if (capacity < 0.0) {
        throw invalid_argument("capacity must be nonnegative");
    }
    if (flow < 0.0 || flow > capacity) {
        throw invalid_argument("flow must be nonnegative and less than or equal to capacity");
    }
    _v = v;
    _w = w;
    _capacity = capacity;
    _flow = flow;
}

FlowEdge::FlowEdge(FlowEdge &e) {
    _v = e._v;
    _w = e._w;
    _capacity = e._capacity;
    _flow = e._flow;
}

FlowEdge::FlowEdge(FlowEdge &&e) {
    _v = e._v;
    _w = e._w;
    _capacity = e._capacity;
    _flow = e._flow;
}

FlowEdge &FlowEdge::operator=(FlowEdge &&e) {
    if (this != &e) {
        _v = e._v;
        _w = e._w;
        _capacity = e._capacity;
        _flow = e._flow;
    }
    return *this;
}

FlowEdge &FlowEdge::operator=(FlowEdge &e) {
    if (this != &e) {
        _v = e._v;
        _w = e._w;
        _capacity = e._capacity;
        _flow = e._flow;
    }
    return *this;
}

int FlowEdge::from() const {
    return _v;
}

int FlowEdge::to() const {
    return _w;
}

double FlowEdge::capacity() const {
    return _capacity;
}

double FlowEdge::flow() const {
    return _flow;
}

int FlowEdge::other(int vertex) const {
    if (vertex == _v) return _w;
    else if (vertex == _w) return _v;
    else throw invalid_argument("invalid vertex");
}

double FlowEdge::residualCapacityTo(int vertex) const {
    if (vertex == _v) return _flow;
    else if (vertex == _w) return _capacity - _flow;
    else throw invalid_argument("invalid vertex");
}

void FlowEdge::addResidualFlowTo(int v, double delta) {
    if (delta < 0.0) {
        throw invalid_argument("delta must be nonnegative");
    }
    if (v == _v) _flow -= delta;
    else if (v == _w) _flow += delta;
    else throw invalid_argument("invalid vertex");
    if (abs(_flow) <= FLOATING_POINT_EPSILON) _flow = 0.0;
    if (abs(_flow - _capacity) <= FLOATING_POINT_EPSILON) _flow = _capacity;
    if (_flow < 0.0) {
        throw invalid_argument("flow is negative");
    }
    if (_flow > _capacity) {
        throw invalid_argument("flow exceeds capacity");
    }
}

string FlowEdge::toString() const {
    stringstream ss;
    ss << _v << "->" << _w << " " << _flow << "/" << _capacity;
    return ss.str();
}

ostream &operator<<(ostream &os, const FlowEdge &e) {
    os << e.toString();
    return os;
}

class CapacityComparator {
public:
    bool operator()(const FlowEdge& e1, const FlowEdge& e2) const {
        return e1.capacity() < e2.capacity();
    }
};
class FlowComparator {
public:
    bool operator()(const FlowEdge& e1, const FlowEdge& e2) const {
        return e1.flow() < e2.flow();
    }
};



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
        int _V;
        int _E;
        int _source;
        int _sink;
    public:
        Bag<FlowEdge*> *_adj;
        FlowNetwork(int V);
        FlowNetwork(int V, int E, int source, int sink);
        FlowNetwork(istream &in);
        ~FlowNetwork();
        void validateVertex(int v) const;
        void addEdge(FlowEdge* e);
        void addEdge(int v, int w, double capacity);
        void removeEdge(FlowEdge *pEdge);
        bool containsEdge(FlowEdge *pEdge);
        int V() const;
        int E() const;
        int source() const { return _source; }
        int sink() const { return _sink; }
        void setSource(int source) { _source = source; }
        void setSink(int sink) { _sink = sink; }
        Bag<FlowEdge*>::Iterator adj(int v);
        Bag<FlowEdge*>::Iterator edges();

        string toString() const;
        friend ostream &operator<<(ostream &os, const FlowNetwork &G);

        // inner class for iterating over the edges
        class FlowNetworkIterator : public Bag<FlowEdge>::Iterator {
        public:
            FlowNetworkIterator(Bag<FlowEdge>::Iterator *it) : Bag<FlowEdge>::Iterator(
                    reinterpret_cast<Node<struct FlowEdge> *>(it)) {}
            FlowNetworkIterator(const FlowNetworkIterator &it) : Bag<FlowEdge>::Iterator(it) {}
            FlowNetworkIterator &operator=(const FlowNetworkIterator &it) {
                Bag<FlowEdge>::Iterator::operator=(it);
                return *this;
            }
            FlowNetworkIterator &operator++() {
                Bag<FlowEdge>::Iterator::operator++();
                return *this;
            }
            FlowNetworkIterator operator++(int) {
                FlowNetworkIterator it = *this;
                operator++();
                return it;
            }
        };



    };


FlowNetwork::FlowNetwork(int V) {
    if (V < 0) throw "Number of vertices must be nonnegative";
    _V = V;
    _E = 0;
    _adj = new Bag<FlowEdge*>[V];
    for (int v = 0; v < V; v++) {
        _adj[v] = Bag<FlowEdge*>();
    }
}

FlowNetwork::FlowNetwork(int V, int E, int source, int sink) {
    if (V < 0) throw "Number of vertices must be nonnegative";
    if (E < 0) throw "Number of edges must be nonnegative";
    if (source < 0 || source >= V) throw "Invalid source";
    if (sink < 0 || sink >= V) throw "Invalid sink";
    _V = V;
    _E = 0;
    _source = source;
    _sink = sink;
    _adj = new Bag<FlowEdge*>[V];
    for (int v = 0; v < V; v++) {
        _adj[v] = Bag<FlowEdge*>();
    }
    // initialize a random flow network with V vertices and E edges. The
    // capacities are between 0 and 99 and the flow values are 0;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> vert(0, V - 1);
    std::uniform_int_distribution<> dis(0, 99); // for capacities
    for (int i = 0; i < E; i++) {
        // randomly generate the vertices from 0 to V-1
        int v = vert(gen);
        int w = vert(gen);
        int c = dis(gen);
        if (v == w) { // no self loops
            i--;
            continue;
        }
        while (v == _sink) v = vert(gen); // no output edges from the sink
        while (w == _source) w = vert(gen); // no input edges to the source
        FlowEdge *e = new FlowEdge(v, w, c);
        if (containsEdge(e)) {
            delete e;
            i--;
            continue;
        } else {
            addEdge(e);
        }
    }
}

FlowNetwork::FlowNetwork(istream &in) {
    in >> _V;
    if (_V < 0) throw "Number of vertices must be nonnegative";
    in >> _E;
    if (_E < 0) throw "Number of edges must be nonnegative";
    _adj = new Bag<FlowEdge*>[V()];
    for (int v = 0; v < V(); v++) {
        _adj[v] = Bag<FlowEdge*>();
    }
    for (int i = 0; i < E(); i++) {
        int v, w;
        double c;
        in >> v >> w >> c;
        if (v < 0 || v >= V()) throw "vertex " + to_string(v) + " is not between 0 and " + to_string(V() - 1);
        if (w < 0 || w >= V()) throw "vertex " + to_string(w) + " is not between 0 and " + to_string(V() - 1);
        addEdge(new FlowEdge(v, w, c));
    }
}

FlowNetwork::~FlowNetwork() {
    unordered_set<FlowEdge*> edges;
    // before deleting the edge add to a set to avoid duplicate deletion
    for (int v = 0; v < V(); v++) {
        for (int i = 0; i < _adj[v].size(); i++) {
            edges.insert(_adj[v].get(i));
        }
    }
    for (auto e : edges) {
        if (e != nullptr) delete e;
    }
    delete[] _adj;
}

void FlowNetwork::validateVertex(int v) const {
    if (v < 0 || v >= V()) throw "vertex " + to_string(v) +
                                 " is not between 0 and " + to_string(V() - 1);
}
void FlowNetwork::addEdge(FlowEdge* e) {
    int v = e->from();
    int w = e->to();
    validateVertex(v);
    validateVertex(w);
    _adj[v].add(e);
    _adj[w].add(e);
    _E++;
}

int FlowNetwork::V() const {
    return _V;
}

int FlowNetwork::E() const {
    return _E;
}

Bag<FlowEdge*>::Iterator FlowNetwork::adj(int v)  {
    validateVertex(v);
    return _adj[v].begin();
}

Bag<FlowEdge*>::Iterator FlowNetwork::edges()  {
    Bag<FlowEdge*> bag;
    for (int v = 0; v < V(); v++) {
        for (FlowEdge* e : _adj[v]) {
            if (e->from() != v) bag.add(e);
        }
    }
    return bag.begin();
}

string FlowNetwork::toString() const {
    stringstream ss;
    ss << V() << " " << E() << endl;
    for (int v = 0; v < V(); v++) {
        ss << "  " << v << ": ";
        for (Bag<FlowEdge*>::Iterator it = _adj[v].begin(); it != _adj[v].end(); ++it) {
            if ((*it)->to() != v) ss << *(*it) << " ";
        }
        ss << endl;
    }
    return ss.str();
}

ostream &operator<<(ostream &os, const FlowNetwork &G) {
    os << G.toString();
    return os;
}

void FlowNetwork::removeEdge(FlowEdge *pEdge) {
    // remove all occurrences of pEdge from _adj
    int v = pEdge->from();
    int w = pEdge->to();
    _adj[v].remove(pEdge);
    _adj[w].remove(pEdge);
    delete pEdge;
    _E--;
}

bool FlowNetwork::containsEdge(FlowEdge *pEdge) {
    // remove all occurrences of pEdge from _adj
    int v = pEdge->from();
    int w = pEdge->to();
    for (Bag<FlowEdge*>::Iterator it = _adj[v].begin(); it != _adj[v].end(); it++) {
        // check that the to and from vertices are the same
        if ((*it)->to() == w && (*it)->from() == v) {
            return true;
        }
    }
    return false;
}

void FlowNetwork::addEdge(int v, int w, double capacity) {
    if (v < 0 || v >= V()) throw "vertex " + to_string(v) + " is not between 0 and " + to_string(V() - 1);
    if (w < 0 || w >= V()) throw "vertex " + to_string(w) + " is not between 0 and " + to_string(V() - 1);
    FlowEdge *e = new FlowEdge(v, w, capacity);
    if (containsEdge(e)) {
        delete e;
        return;
    } else {
        addEdge(e);
    }
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
    int  _V; // number of vertices
    int _source; // source vertex
    int _sink; // sink vertex
    // marked[v] = true iff s->v path in residual graph
    vector<bool> marked;
    // edgeTo[v] = last edge on shortest s->v path
    vector<FlowEdge*> edgeTo;
    double _value;
    FlowNetwork* _network;

public:
    FordFulkerson(FlowNetwork &G, int s, int t);
    double value();
    bool inCut(int v);
    void validate(int v);
    bool hasAugmentingPath(FlowNetwork& G, int s, int t);
    double excess(FlowNetwork& G, int v);
    bool isFeasible(FlowNetwork& G, int s, int t);
    bool check(FlowNetwork& G, int s, int t);
    string toString();
    // overload << operator
    friend ostream& operator<<(ostream& os, FordFulkerson& ff);
    ~FordFulkerson() {
        _network = nullptr;
        delete _network;
    }
};

FordFulkerson::FordFulkerson(FlowNetwork &G, int s, int t) {
    _V = G.V();
    validate(s);
    validate(t);
    if (s == t) throw runtime_error("s == t");
    _network = &G;
    _source = s;
    _sink = t;
    _value = 0.0;
    if (!isFeasible(G, s, t)) throw runtime_error("Initial flow is infeasible");
    edgeTo = vector<FlowEdge*>(_V);
    marked = vector<bool>(_V);
    for (int v = 0; v < _V; v++) {
        edgeTo[v] = nullptr;
        marked[v] = false;
    }
    _value= excess(G, t);
    // compute maximum flow and minimum cut
    while (hasAugmentingPath(G, s, t)) {
        double bottle = std::numeric_limits<double>::max();
        for (int v = t; v != s; v = edgeTo[v]->other(v)) {
            bottle = min(bottle, edgeTo[v]->residualCapacityTo(v));
        }
        // augment flow
        for (int v = t; v != s; v = edgeTo[v]->other(v)) {
            edgeTo[v]->addResidualFlowTo(v, bottle);
        }
        _value += bottle;
    }

    // check optimality conditions
    if (!isFeasible(G, s, t)) {
        cout << "*** Flow is infeasible ***" << endl;
    } else if (!check(G, s, t)) {
        cout << "*** check() detects a problem ***" << endl;
    } else {
        cout << "*** Flow is optimal ***" << endl;
    }

}

double FordFulkerson::value() {
    return _value;
}

bool FordFulkerson::inCut(int v) {
    validate(v);
    return marked[v];
}

void FordFulkerson::validate(int v) {
    if (v < 0 || v >= _V)
        throw "vertex " + to_string(v) + " is not between 0 and " +
              to_string(_V - 1);
}
bool FordFulkerson::hasAugmentingPath(FlowNetwork &G, int s, int t) {
    edgeTo = vector<FlowEdge*>(_V);
    marked = vector<bool>(_V);
    for (int v = 0; v < _V; v++) {
        edgeTo[v] = nullptr;
        marked[v] = false;
    }
    queue<int> q;
    q.push(s);
    marked[s] = true;
    while (!q.empty() && !marked[t]) {
        int v = q.front();
        q.pop();
        for (Bag<FlowEdge*>::Iterator e = G._adj[v].begin(); e !=
                                                             G._adj[v].end(); ++e) {
            int w = (*e)->other(v);
            if ((*e)->residualCapacityTo(w) > 0) {
                if (!marked[w]) {
                    edgeTo[w] = *e;
                    marked[w] = true;
                    this->marked[w] = true;
                    q.push(w);
                }
            }
        }
    }
    return marked[t];
}

double FordFulkerson::excess(FlowNetwork &G, int v) {
    double excess = 0.0;
    for (Bag<FlowEdge*>::Iterator e = G._adj[v].begin(); e != G._adj[v].end(); ++e) {
        if ((*e)->from() == v) {
            excess -= (*e)->flow();
        } else {
            excess += (*e)->flow();
        }
    }
    return excess;
}

bool FordFulkerson::isFeasible(FlowNetwork &G, int s, int t) {
    for (int v = 0; v < G.V(); v++) {
        for (Bag<FlowEdge*>::Iterator e = G._adj[v].begin(); e != G._adj[v].end(); ++e) {
            if ((*e)->flow() < -FLOATING_POINT_EPSILON ||
                (*e)->flow() > (*e)->capacity() + FLOATING_POINT_EPSILON) {
                std::cout << "Edge does not satisfy capacity constraints: " <<
                          *e << std::endl;
                return false;
            }
        }
    }
    if (abs(_value + excess(G, s)) > FLOATING_POINT_EPSILON) {
        std::cout << "Excess at source: " << s << " = " << excess(G, s) <<
                  ", Max Flow = " << _value << std::endl;
        return false;
    }
    if (abs(_value - excess(G, t)) > FLOATING_POINT_EPSILON) {
        std::cout << "Excess at sink: " << t << " = " << excess(G, t) <<
                  ", Max Flow = " << _value << std::endl;
        return false;
    }
    for (int v = 0; v < G.V(); v++) {
        if (v == s || v == t) continue;
        else if (abs(excess(G, v)) > FLOATING_POINT_EPSILON) {
            std::cout << "Net flow out of " << v <<
                      " doesn't equal zero" << std::endl;
            return false;
        }
    }
    return true;
}

bool FordFulkerson::check(FlowNetwork &G, int s, int t) {
    if (!isFeasible(G, s, t)) {
        std::cout << "Flow is infeasible" << std::endl; return false;
    }
    if (!inCut(s)) {
        std::cout << "source " << s << " is not on source side of min cut" <<
                  std::endl; return false;
    }
    if (inCut(t)) {
        std::cout << "sink " << t << " is on source side of min cut" <<
                  std::endl; return false;
    }
    double mincutValue = 0.0;
    for (int v = 0; v < G.V(); v++) {
        for (Bag<FlowEdge*>::Iterator e = G._adj[v].begin(); e != G._adj[v].end(); ++e) {
            if ((v == (*e)->from() && inCut((*e)->from()) && !inCut((*e)->to()))) {
                mincutValue += (*e)->capacity();
            }
        }
    }
    if (abs(mincutValue - _value) > FLOATING_POINT_EPSILON) {
        std::cout << "Max flow value = " << _value << ", min cut value = "
                  << mincutValue << std::endl; return false;
    }
    return true;
}

string FordFulkerson::toString() {
    string s = "";
    s += "Max flow from  " + to_string(_source) + " to " + to_string(_sink) + "\n";
    for (int v = 0; v < _V; v++) {
        for (Bag<FlowEdge*>::Iterator it = _network->_adj[v].begin(); it !=
                                                                      _network->_adj[v].end(); ++it) {
            if ((*it)->from() == v&& (*it)->flow() > 0) {
                s+= " " + (*it)->toString() + "\n";
            }
        }
    }
    s += "Min cut: ";
    for (int v = 0; v < _V; v++) {
        if (inCut(v)) s += to_string(v) + " ";
    }
    s += "\nMax Flow Value: " + to_string(_value);
    s += "\n";
    return s;
}

ostream &operator<<(ostream &os, FordFulkerson &ff) {
    os << ff.toString();
    return os;
}

#endif //FLOW_NETWORK_FORDFULKERSON_H