//
// Created by Ryan.Zurrin001 on 1/13/2022.
//

#ifndef GRAPH_CPP_FORDFULKERSON_H
#define GRAPH_CPP_FORDFULKERSON_H
#include "FlowNetwork.h"

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
#endif //GRAPH_CPP_FORDFULKERSON_H
