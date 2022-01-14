//
// Created by Ryan.Zurrin001 on 1/13/2022.
//

#ifndef GRAPH_CPP_FORD_FULKERSON_H
#define GRAPH_CPP_FORD_FULKERSON_H
#include "Flow_Network.h"
#define DEBUG 1

class Ford_Fulkerson {

    Flow_Network & network_;
    int V_;
    int value_;
    vector<bool> visited_;
    vector<int> path_;
    vector<Flow_Network> edges_;

public:
    Ford_Fulkerson(Flow_Network& network, int s, int t);
    ~Ford_Fulkerson();
    int get_value();
    vector<int> get_path();
    void ford_fulkerson_algorithm();
    void find_path(int v);
    void augment(int v);
    void find_max_flow();
    bool validate(int v);
    bool inCut(int v);
    bool hasAugmentingPath(int v);
    double excess(int v);
    bool isFeasible(int s, int t);
    bool check(int s, int t);
    void print_max_flow();

};

Ford_Fulkerson::Ford_Fulkerson(Flow_Network& network, int s, int t) :
network_(network){
    if (DEBUG) {
        cout << "Ford_Fulkerson Constructor" << endl;
    }
    V_ = network.V();
    validate(s);
    validate(t);
    if (s == t) {
        throw invalid_argument("Source and sink cannot be the same");
    }
    visited_.resize(network.E());
    path_.resize(network.E());
    value_ = excess(t);

    for (int i = 0; i < V_; ++i) {
        edges_.push_back(Flow_Network(network.E()));
    }
    // loop through all edges


    while (hasAugmentingPath(s)) {
        double bottleneck = numeric_limits<double>::max();
        for (int i = 0; i < path_.size(); i++) {
            if (path_[i] != -1) {
                bottleneck = min(bottleneck, network.get_capacity(path_[i]));
        }
        }
        for (int i = 0; i < path_.size(); i++) {
            if (path_[i] != -1) {
                network.set_capacity(path_[i], network.get_capacity(path_[i]) - bottleneck);
            }
        }
        value_ += bottleneck;
    }
}

int Ford_Fulkerson::get_value() {
    if (DEBUG) cout << "Ford_Fulkerson get_value" << endl;
        return value_;
}

vector<int> Ford_Fulkerson::get_path() {
    if (DEBUG) cout << "Ford_Fulkerson::get_path()" << endl;
        return path_;
}

void Ford_Fulkerson::find_max_flow() {
    if (DEBUG) cout << "Ford_Fulkerson::find_max_flow()" << endl;
        while (true) {
        for (int i = 0; i < V_; i++) {
            visited_[i] = false;
            path_[i] = -1;
        }
        find_path(0);
        if (visited_[V_ - 1]) {
            augment(0);
        } else {
            break;
        }
    }
    if (DEBUG) cout << "Ford_Fulkerson::find_max_flow() End" << endl;
}

bool Ford_Fulkerson::validate(int v) {
    if (DEBUG) cout << "Ford_Fulkerson::validate()" << endl;
        if (v < 0 || v >= V_) {
        throw invalid_argument("Vertex is out of bounds");
    }
    if (DEBUG) cout << "Ford_Fulkerson::validate() End" << endl;
    return true;
}

void Ford_Fulkerson::print_max_flow() {
    if (DEBUG) cout << "Ford_Fulkerson::print_max_flow()" << endl;
    cout << "Max flow: " << value_ << endl;
    for (int i = 0; i < V_; i++) {
        for (int j = 0; j < V_; j++) {
            cout << edges_[i].edge(i, j) << " ";
        }
        cout << endl;
    }
    if (DEBUG) cout << "Ford_Fulkerson::print_max_flow() End" << endl;
}

void Ford_Fulkerson::ford_fulkerson_algorithm() {
    cout << "Ford Fulkerson Algorithm" << endl;
    cout << "Max flow: " << value_ << endl;
    for (int i = 0; i < V_; i++) {
        for (int j = 0; j < V_; j++) {
            cout << edges_[i].edge(i, j) << " ";
        }
        cout << endl;
    }
    if (DEBUG) cout << "Ford_Fulkerson::ford_fulkerson_algorithm() End" << endl;
}

bool Ford_Fulkerson::inCut(int v) {
        if (DEBUG) cout << "Ford_Fulkerson::inCut()" << endl;
    validate(v);
    if (DEBUG) cout << "visited_[v]: " << visited_[v] << endl;
        return visited_[v];
}

bool Ford_Fulkerson::hasAugmentingPath(int v) {
    if (DEBUG) cout << "Ford_Fulkerson::hasAugmentingPath()" << endl;
    validate(v);
    for (int i = 0; i < V_; i++) {
        if (edges_[v].edgeCost(v, i) > 0 && !visited_[i]) {
            return true;
        }
    }
    if (DEBUG) cout << "Ford_Fulkerson::hasAugmentingPath() End" << endl;
    return false;
}

double Ford_Fulkerson::excess(int v) {
    if (DEBUG) cout << "Ford_Fulkerson::excess()" << endl;
    validate(v);
    double excess = 0;
    for (int i = 0; i < V_; i++) {
        if (v == edges_[i].get_source()) {
            excess -= edges_[i].edgeCost(v, i);
        } else if (v == edges_[i].get_destination()) {
            excess += edges_[i].edgeCost(v, i);
        }
            excess -= edges_[i].edgeCost(v, i);
        } else if (v == edges_[i].get_destination()) {
            excess += edges_[i].edgeCost(v, i);
        }
    }
}

bool Ford_Fulkerson::isFeasible(int s, int t) {
    if (DEBUG) cout << "Ford_Fulkerson::isFeasible()" << endl;
    validate(s);
    validate(t);
    if (excess(s) < 0) {
        return false;
    }
    for (int i = 0; i < V_; i++) {
        if (i != s && i != t && excess(i) > 0) {
            return false;
        }
    }
    if (DEBUG) cout << "Ford_Fulkerson::isFeasible() End" << endl;
    return true;
}

bool Ford_Fulkerson::check(int s, int t) {
    if (DEBUG) cout << "Ford_Fulkerson::check()" << endl;
    if (!isFeasible(s, t)) {
        cout << "The graph is not feasible" << endl;
        return false;
    }
    for (int i = 0; i < V_; i++) {
        if (i != s && i != t && excess(i) > 0) {
            cout << "The graph is not bipartite" << endl;
            return false;
        }
    }
    if (DEBUG) cout << "Ford_Fulkerson::check() End" << endl;
    return true;
}

void Ford_Fulkerson::augment(int v) {
    if (DEBUG) cout << "Ford_Fulkerson::augment()" << endl;
    visited_[v] = true;
    for (int i = 0; i < V_; i++) {
        if (edges_[v].edgeCost(v, i) > 0 && !visited_[i]) {
            if (excess(i) > 0) {
                edges_[v].addEdge(v, i, edges_[v].edgeCost(v, i));
                edges_[i].addEdge(i, v, edges_[v].edgeCost(v, i));
                edges_[v].addEdge(v, i, 0);
                edges_[i].addEdge(i, v, 0);
                augment(i);
            }
        }
    }
    if (DEBUG) cout << "Ford_Fulkerson::augment() End" << endl;

}

void Ford_Fulkerson::find_path(int v) {
    if (DEBUG) cout << "Ford_Fulkerson::find_path(" << v << ")" << endl;
        visited_[v] = true;
    for (int i = 0; i < V_; i++) {
        if (edges_[v].edgeCost(v, i) > 0 && !visited_[i]) {
            if (excess(i) > 0) {
                path_[i] = v;
                find_path(i);
            }
        }
    }
    if (DEBUG) cout << "Ford_Fulkerson::find_path(" << v << ") End" << endl;
}

Ford_Fulkerson::~Ford_Fulkerson() {
    if (DEBUG) cout << "Ford_Fulkerson::~Ford_Fulkerson()" << endl;
    edges_.clear();
    visited_.clear();
    path_.clear();
    if (DEBUG) cout << "Ford_Fulkerson::~Ford_Fulkerson() End" << endl;
}


#endif //GRAPH_CPP_FORD_FULKERSON_H
