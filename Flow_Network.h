//
// Created by Ryan.Zurrin001 on 1/13/2022.
//

#ifndef GRAPH_CPP_FLOW_NETWORK_H
#define GRAPH_CPP_FLOW_NETWORK_H
#include "Flow_Edge.h"

class Flow_Network {
public:
    Flow_Network(int V);
    Flow_Network(int V, int E);
    Flow_Network(const Flow_Network& other);
    Flow_Network& operator=(const Flow_Network& other);
    ~Flow_Network();

    int V() const;
    int E() const;
    int flow() const;
    void addEdge(Flow_Edge* e);
    void addEdge(int v, int w, double capacity);
    void addEdge(int v, int w, double capacity, double cost);
    bool validEdge(int v, int w) const;
    // return vector iterator to the edges
    vector<Flow_Edge*>::iterator adj(int v);
    vector<Flow_Edge*>::iterator edges() const;
    // get the edge
    Flow_Edge* edge(int v, int w) const;
    double edgeCost(int v, int w) const;
    // get val at the edge
    // overload the << operator to print the graph
    friend ostream& operator<<(ostream& os, const Flow_Network& g);
    // overload the subscript operator to return the edge
    Flow_Edge& operator[](int v);
    // overload the () operator to return the edge
    Flow_Edge& operator()(int v);
    Flow_Edge& operator()(int v, int w);

    Flow_Edge *edge(int i);
    vector<Flow_Edge*>* _adj;


   double get_capacity(int &i);

    void set_capacity(int &i, double d);

private:
    int _V;
    int _E;


};

Flow_Network::Flow_Network(int V) {
    _V = V;
    _E = 0;
    _adj = new vector<Flow_Edge*>[_V];
}

Flow_Network::Flow_Network(int V, int E) {
    if (V < 0) {
        throw std::invalid_argument("Number of vertices must be nonnegative");
    }
    if (E < 0) {
        throw std::invalid_argument("Number of edges must be nonnegative");
    }
    _V = V;
    _E = 0;
    _adj = new vector<Flow_Edge*>[_V];
    // initialize a random flow network with V vertices and E edges. The capacities
    // are integer values between 0 and 99 and the flow values are all 0.
    srand(time(NULL));
    for (int i = 0; i < E; i++) {
        int v = rand() % V;
        int w = rand() % V;
        int c = rand() % 100;
        while (v == w) {
            w = rand() % V;
        }
        // check to see if the edge already exists
        if (validEdge(v, w)) {
            i--;
            continue;
        }
        addEdge(v, w, c);
    }
}

Flow_Network::Flow_Network(const Flow_Network &other) {
    _V = other._V;
    _E = other._E;
    _adj = new vector<Flow_Edge*>[_V];
    for (int v = 0; v < _V; v++) {
        for (Flow_Edge* e : other._adj[v]) {
            _adj[v].push_back(new Flow_Edge(*e));
        }
    }
}

Flow_Network &Flow_Network::operator=(const Flow_Network &other) {
    if (this != &other) {
        _V = other._V;
        _E = other._E;
        _adj = new vector<Flow_Edge*>[_V];
        for (int v = 0; v < _V; v++) {
            for (Flow_Edge* e : other._adj[v]) {
                _adj[v].push_back(new Flow_Edge(*e));
            }
        }
    }
    return *this;
}

Flow_Network::~Flow_Network() {
    for (int v = 0; v < _V; v++) {
        for (Flow_Edge* e : _adj[v]) {
            delete e;
        }
    }
    delete[] _adj;
}

int Flow_Network::V() const {
    return _V;
}

int Flow_Network::flow() const {
    int flow = 0;
    for (int v = 0; v < _V; v++) {
        for (Flow_Edge* e : _adj[v]) {
            flow += e->flow();
        }
    }
    return flow;
}

int Flow_Network::E() const {
    return _E;
}

void Flow_Network::addEdge(Flow_Edge *e) {
    int v = e->from();
    int w = e->to();
    if (v < 0 || v >= _V) {
        cout << "vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1);
        return;
    }
    if (w < 0 || w >= _V) {
        cout << "vertex " + std::to_string(w) + " is not between 0 and " + std::to_string(_V - 1) << endl;
        return;
    }
    if (v == w) {
        cout << "vertex " + std::to_string(v) + " is equal to vertex " + std::to_string(w) + "\n";
        return;
    }
    if (validEdge(v, w)) {
        cout << "edge (" + std::to_string(v) + ", " + std::to_string(w) + ") is a duplicate" << endl;
        return;
    }
    _adj[v].push_back(e);
    _E++;
}

void Flow_Network::addEdge(int v, int w, double capacity) {
    if (v < 0 || v >= _V) {
        cout << "vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1);
        return;
    }
    if (w < 0 || w >= _V) {
        cout << "vertex " + std::to_string(w) + " is not between 0 and " + std::to_string(_V - 1) << endl;
        return;
    }
    if (v == w) {
        cout << "vertex " + std::to_string(v) + " is equal to vertex " + std::to_string(w) + "\n";
        return;
    }
    if (validEdge(v, w)) {
        cout << "edge (" + std::to_string(v) + ", " + std::to_string(w) + ") is a duplicate" << endl;
        return;
    }
    _adj[v].push_back(new Flow_Edge(v, w, capacity));
    _E++;
}

void Flow_Network::addEdge(int v, int w, double capacity, double cost) {
    if (v < 0 || v >= _V) {
        cout << "vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1);
        return;
    }
    if (w < 0 || w >= _V) {
        cout << "vertex " + std::to_string(w) + " is not between 0 and " + std::to_string(_V - 1) << endl;
        return;
    }
    if (v == w) {
        cout << "vertex " + std::to_string(v) + " is equal to vertex " + std::to_string(w) + "\n";
        return;
    }
    if (validEdge(v, w)) {
        cout << "edge (" + std::to_string(v) + ", " + std::to_string(w) + ") is a duplicate" << endl;
        return;
    }
    _adj[v].push_back(new Flow_Edge(v, w, capacity, cost));
    _E++;
}

vector<Flow_Edge *>::iterator Flow_Network::adj(int v) {
    if (v < 0 || v >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1));
    }
    return _adj[v].begin();
}

bool Flow_Network::validEdge(int v, int w) const {
    for (Flow_Edge* e : _adj[v]) {
        if (e->to() == w) {
            return true;
        }
    }
    return false;
}

vector<Flow_Edge *>::iterator Flow_Network::edges() const {
    vector<Flow_Edge *> edges;
    for (int v = 0; v < _V; v++) {
        for (Flow_Edge* e : _adj[v]) {
            edges.push_back(e);
        }
    }
    return edges.begin();
}

ostream &operator<<(ostream &os, const Flow_Network &g) {
    os << "Flow_Network(" << g.V() << " vertices, " << g.E() << " edges)" << endl;
    for (int v = 0; v < g.V(); v++) {
        os << "  " << v << ": ";
        for (Flow_Edge* e : g._adj[v]) {
            os << *e << "  ";
        }
        os << endl;
    }
    return os;
}

Flow_Edge &Flow_Network::operator[](int v) {
    if (v < 0 || v >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1));
    }
    return *_adj[v][0];
}

Flow_Edge &Flow_Network::operator()(int v) {
    if (v < 0 || v >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1));
    }
    return *_adj[v][0];
}

Flow_Edge &Flow_Network::operator()(int v, int w) {
    if (v < 0 || v >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1));
    }
    if (w < 0 || w >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(w) + " is not between 0 and " + std::to_string(_V - 1));
    }
    for (Flow_Edge* e : _adj[v]) {
        if (e->to() == w) {
            return *e;
        }
    }
    throw std::invalid_argument("edge (" + std::to_string(v) + ", " + std::to_string(w) + ") is not in the graph");
}

Flow_Edge *Flow_Network::edge(int v, int w) const {
    if (v < 0 || v >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1));
    }
    if (w < 0 || w >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(w) + " is not between 0 and " + std::to_string(_V - 1));
    }
    for (Flow_Edge* e : _adj[v]) {
        if (e->to() == w) {
            return e;
        }
    }
    return nullptr;
}


Flow_Edge *Flow_Network::edge(int i) {
    if (i < 0 || i >= _E) {
        throw std::invalid_argument("edge " + std::to_string(i) + " is not between 0 and " + std::to_string(_E - 1));
    }
    int count = 0;
    for (int v = 0; v < _V; v++) {
        for (Flow_Edge* e : _adj[v]) {
            if (count == i) {
                return e;
            }
            count++;
        }
    }
    throw std::invalid_argument("edge " + std::to_string(i) + " is not in the graph");
}

double Flow_Network::edgeCost(int v, int w) const {
    if (v < 0 || v >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(v) + " is not between 0 and " + std::to_string(_V - 1));
    }
    if (w < 0 || w >= _V) {
        throw std::invalid_argument("vertex " + std::to_string(w) + " is not between 0 and " + std::to_string(_V - 1));
    }
    for (Flow_Edge* e : _adj[v]) {
        if (e->to() == w) {
            return e->capacity() - e->flow();
        }
    }
    throw std::invalid_argument("edge (" + std::to_string(v) + ", " + std::to_string(w) + ") is not in the graph");
}


double  Flow_Network::get_capacity(int &i) {
    if (i < 0 || i >= _E) {
        throw std::invalid_argument("edge " + std::to_string(i) + " is not between 0 and " + std::to_string(_E - 1));
    }
    int count = 0;
    for (int v = 0; v < _V; v++) {
        for (Flow_Edge* e : _adj[v]) {
            if (count == i) {
                return e->residualCapacityTo(v);
            }
            count++;
        }
    }
    throw std::invalid_argument("edge " + std::to_string(i) + " is not in the graph");
}

void Flow_Network::set_capacity(int &i, double d) {
    // adds resodiial flow to the edge
    if (i < 0 || i >= _E) {
        throw std::invalid_argument("edge " + std::to_string(i) + " is not between 0 and " + std::to_string(_E - 1));
    }
    int count = 0;
    for (int v = 0; v < _V; v++) {
        for (Flow_Edge* e : _adj[v]) {
            if (count == i) {
                e->addResidualFlowTo(v, d);
                return;
            }
            count++;
        }
    }
    throw std::invalid_argument("edge " + std::to_string(i) + " is not in the graph");
}







//
//class Flow_Network {
//public:
//    std::vector<std::shared_ptr<Flow_Edge>>* _adj;
//    Flow_Network(int V);
//    Flow_Network(int V, int E);
//    // move semantics
//    Flow_Network(Flow_Network &&other);
//    Flow_Network &operator=(vector<shared_ptr<Flow_Edge>> other);
//    // copy semantics
//    Flow_Network(const Flow_Network &other);
//    Flow_Network &operator=(const Flow_Network &other);
//    // constructor for taking in input from input stream
//    Flow_Network(std::istream& in);
//    ~Flow_Network();
//
//    int V() const;
//    int E() const;
//    void validateVertex(int v);
//    void addEdge(std::shared_ptr<Flow_Edge> e);
//    void addEdge(int v, int w, double capacity);
//    void addEdge(Flow_Edge e);
//    // an iterator for adj
//    std::vector<std::shared_ptr<Flow_Edge>>::iterator adj(int v);
//    std::shared_ptr<Flow_Edge> edge(int v, int w);
//    int getEdge(int v, int w);
//    std::vector<std::shared_ptr<Flow_Edge>>::iterator edges();
//    string toString()const;
//    // overload << operator
//    friend std::ostream& operator<<(std::ostream& out, const Flow_Network& G);
//    // implement an iterator for edges
//    class iterator {
//    public:
//        iterator(std::vector<std::shared_ptr<Flow_Edge>>::iterator it) {
//            _it = it;
//        }
//        iterator(const iterator &other) {
//            _it = other._it;
//        }
//        iterator &operator=(const iterator &other) {
//            _it = other._it;
//            return *this;
//        }
//        iterator &operator++() {
//            ++_it;
//            return *this;
//        }
//        iterator operator++(int) {
//            iterator temp = *this;
//            ++_it;
//            return temp;
//        }
//        bool operator==(const iterator &other) const {
//            return _it == other._it;
//        }
//        bool operator!=(const iterator &other) {
//            return _it != other._it;
//        }
//        std::shared_ptr<Flow_Edge> operator*() {
//            return *_it;
//        }
//        std::shared_ptr<Flow_Edge> *operator->() {
//            return &(*_it);
//        }
//    private:
//        std::vector<std::shared_ptr<Flow_Edge>>::iterator _it;
//    };
//    iterator begin() {
//        return iterator(_adj->begin());
//    }
//    iterator end() {
//        return iterator(_adj->end());
//    }
//
//private:
//    int _V;
//    int _E;
//};
//
//Flow_Network::Flow_Network(int V) {
//    _V = V;
//    _E = 0;
//    _adj = new std::vector<std::shared_ptr<Flow_Edge>>[_V];
//    for (int i = 0; i < _V; i++) {
//        _adj[i] = std::vector<std::shared_ptr<Flow_Edge>>();
//    }
//}
//
//Flow_Network::Flow_Network(int V, int E) {
//    if (V < 0) {
//        throw std::invalid_argument("Number of vertices must be nonnegative");
//    }
//    if (E < 0) {
//        throw std::invalid_argument("Number of edges must be nonnegative");
//    }
//    _V = V;
//    _E = 0;
//    _adj = new std::vector<std::shared_ptr<Flow_Edge>>[_V];
//    // initialize a random flow network with V vertices and E edges. The capacities
//    // are integer values between 0 and 99 and the flow values are all 0.
//    srand(time(NULL));
//    for (int i = 0; i < E; i++) {
//        int v = rand() % V;
//        int w = rand() % V;
//        int c = rand() % 100;
//        while (v == w) {
//            w = rand() % V;
//        }
//        std::shared_ptr<Flow_Edge> e = std::make_unique<Flow_Edge>(v, w, c);
//        addEdge(std::move(e));
//    }
//}
//
//Flow_Network::Flow_Network(istream &in) {
//    in >> _V >> _E;
//    _adj = new std::vector<std::shared_ptr<Flow_Edge>>[_V];
//    for (int i = 0; i < _E; i++) {
//        int v, w, capacity;
//        in >> v >> w >> capacity;
//        validateVertex(v);
//        validateVertex(w);
//        if (capacity < 0) {
//            throw std::invalid_argument("Edge capacity must be nonnegative");
//        }
//        std::shared_ptr<Flow_Edge> e = std::make_unique<Flow_Edge>(v, w, capacity);
//    }
//}
//
//Flow_Network::~Flow_Network() {
//    delete[] _adj;
//}
//
//int Flow_Network::V() const {
//    return _V;
//}
//
//int Flow_Network::E() const {
//    return _E;
//}
//
//void Flow_Network::validateVertex(int v) {
//    if (v < 0 || v >= _V) {
//        throw std::invalid_argument("vertex " + std::to_string(v) +
//        " is not between 0 and " + std::to_string(_V - 1));
//    }
//}
//
//void Flow_Network::addEdge(std::shared_ptr<Flow_Edge> e) {
//    int v = e->from();
//    int w = e->to();
//    validateVertex(v);
//    validateVertex(w);
//    _adj[v].push_back(e);
//    _adj[w].push_back(e);
//    _E++;
//}
//
//std::vector<std::shared_ptr<Flow_Edge>>::iterator Flow_Network::adj(int v) {
//    validateVertex(v);
//    return _adj[v].begin();
//}
//
//std::shared_ptr<Flow_Edge> Flow_Network::edge(int v, int w) {
//    validateVertex(v);
//    validateVertex(w);
//    for (auto e : _adj[v]) {
//        if (e->to() == w) {
//            return e;
//        }
//    }
//    return NULL;
//}
//
//std::vector<std::shared_ptr<Flow_Edge>>::iterator Flow_Network::edges() {
//    return _adj[0].begin();
//}
//string Flow_Network::toString() const {
//    stringstream ss = stringstream();
//    ss << "V = " << _V << ", E = " << _E << endl;
//    for (int v = 0; v < _V; v++) {
//        ss << v << ": ";
//        for(std::shared_ptr<Flow_Edge> e : _adj[v]) {
//            if (e->to() != v) {
//                ss << e->toString() << ",  ";
//            }
//        }
//        ss << endl;
//    }
//    return ss.str();
//}
//
//std::ostream &operator<<(ostream &out, const Flow_Network &G) {
//    out << G.toString();
//    return out;
//}
//
//Flow_Network::Flow_Network(Flow_Network &&other) {
//    _V = other._V;
//    _E = other._E;
//    _adj = other._adj;
//    other._V = 0;
//    other._E = 0;
//    other._adj = NULL;
//}
//
//Flow_Network &Flow_Network::operator=(vector<shared_ptr<Flow_Edge>> other) {
//    _V = other.size();
//    _E = 0;
//    _adj = new std::vector<std::shared_ptr<Flow_Edge>>[_V];
//    for (int i = 0; i < _V; i++) {
//        _adj[i] = std::vector<std::shared_ptr<Flow_Edge>>();
//    }
//    for (auto e : other) {
//        addEdge(e);
//    }
//    return *this;
//}
//
//Flow_Network::Flow_Network(const Flow_Network &other) {
//    _V = other._V;
//    _E = other._E;
//    _adj = new std::vector<std::shared_ptr<Flow_Edge>>[_V];
//    for (int i = 0; i < _V; i++) {
//        _adj[i] = std::vector<std::shared_ptr<Flow_Edge>>();
//    }
//    for (int v = 0; v < _V; v++) {
//        for (auto e : other._adj[v]) {
//            int w = e->to();
//            int c = e->capacity();
//            std::shared_ptr<Flow_Edge> e2 = std::make_unique<Flow_Edge>(v, w, c);
//            _adj[v].push_back(e2);
//            _adj[w].push_back(e2);
//        }
//    }
//}
//
//Flow_Network &Flow_Network::operator=(const Flow_Network &other) {
//    if (this != &other) {
//        _V = other._V;
//        _E = other._E;
//        delete[] _adj;
//        _adj = new std::vector<std::shared_ptr<Flow_Edge>>[_V];
//        for (int i = 0; i < _V; i++) {
//            _adj[i] = std::vector<std::shared_ptr<Flow_Edge>>();
//        }
//        for (int v = 0; v < _V; v++) {
//            for (auto e : other._adj[v]) {
//                int w = e->to();
//                int c = e->capacity();
//                std::shared_ptr<Flow_Edge> e2 = std::make_unique<Flow_Edge>(v, w, c);
//                _adj[v].push_back(e2);
//                _adj[w].push_back(e2);
//            }
//        }
//    }
//    return *this;
//}
//
//void Flow_Network::addEdge(Flow_Edge e) {
//    int v = e.from();
//    int w = e.to();
//    validateVertex(v);
//    validateVertex(w);
//    _adj[v].push_back(std::make_unique<Flow_Edge>(e));
//    _adj[w].push_back(std::make_unique<Flow_Edge>(e));
//    _E++;
//}
//
//void Flow_Network::addEdge(int v, int w, double capacity) {
//    validateVertex(v);
//    validateVertex(w);
//    if (capacity < 0) {
//        throw std::invalid_argument("Edge capacity must be nonnegative");
//    }
//    _adj[v].push_back(std::make_unique<Flow_Edge>(v, w, capacity));
//    _adj[w].push_back(std::make_unique<Flow_Edge>(v, w, capacity));
//    _E++;
//}
//
//int Flow_Network::getEdge(int v, int w) {
//    validateVertex(v);
//    validateVertex(w);
//    for (auto e : _adj[v]) {
//        if (e->to() == w) {
//            return e->capacity();
//        }
//    }
//    return 0;
//}


#endif //GRAPH_CPP_FLOW_NETWORK_H
