//
// Created by Ryan.Zurrin001 on 1/13/2022.
//

#ifndef GRAPH_CPP_FLOW_EDGE_H
#define GRAPH_CPP_FLOW_EDGE_H
#include <bits/stdc++.h>
using namespace std;
constexpr double FLOATING_POINT_EPSILON = 1e-10;


class Flow_Edge {
public:
    Flow_Edge(int v, int w, double capacity);
    Flow_Edge(int v, int w, double capacity, double flow);
    Flow_Edge(const Flow_Edge& e);
    Flow_Edge& operator=(const Flow_Edge& e);
    Flow_Edge& operator=(Flow_Edge&& e);
    ~Flow_Edge() = default;
    int from() const;
    int to() const;
    int other(int vertex) const;
    double capacity() const;
    double flow() const;
    double residualCapacityTo(int vertex) const;
    void addResidualFlowTo(int vertex, double delta);
    string toString() const;
    // overload the << operator
    friend ostream& operator<<(ostream& out, const Flow_Edge& e);
    friend bool operator==(const Flow_Edge& e1, const Flow_Edge& e2);
    friend bool operator!=(const Flow_Edge& e1, const Flow_Edge& e2);
    friend bool operator<(const Flow_Edge& e1, const Flow_Edge& e2);
    friend bool operator<=(const Flow_Edge& e1, const Flow_Edge& e2);
    friend bool operator>(const Flow_Edge& e1, const Flow_Edge& e2);
    friend bool operator>=(const Flow_Edge& e1, const Flow_Edge& e2);

private:
    int v_;
    int w_;
    double capacity_;
    double flow_;
};
class FlowEdgeCapacityComparator {
public:
    bool operator()(const Flow_Edge &e1, const Flow_Edge &e2) const {
        return e1.capacity() < e2.capacity();
    }
};
class FlowEdgeFlowComparator {
public:
    bool operator()(const Flow_Edge &e1, const Flow_Edge &e2) const {
        return e1.flow() < e2.flow();
    }
};

Flow_Edge::Flow_Edge(int v, int w, double capacity) {
    v_ = v;
    w_ = w;
    capacity_ = capacity;
    flow_ = 0.0;
}

Flow_Edge::Flow_Edge(int v, int w, double capacity, double flow) {
    v_ = v;
    w_ = w;
    capacity_ = capacity;
    flow_ = flow;
}

int Flow_Edge::from() const {
    return v_;
}

int Flow_Edge::to() const {
    return w_;
}


int Flow_Edge::other(int vertex) const {
    if (vertex == v_) return w_;
    else if (vertex == w_) return v_;
    else throw runtime_error("invalid vertex");
}


double Flow_Edge::capacity() const {
    return capacity_;
}

double Flow_Edge::flow() const {
    return flow_;
}

double Flow_Edge::residualCapacityTo(int vertex) const {
    if (vertex == v_) return flow_;
    else if (vertex == w_) return capacity_ - flow_;
    else throw "Illegal endpoint";
}

void Flow_Edge::addResidualFlowTo(int vertex, double delta) {
    if (!(delta > 0)) throw "Delta must be positive";
    if (vertex == v_) flow_ -= delta;
    else if (vertex == w_) flow_ += delta;
    else throw "Illegal endpoint";
    if (abs(flow_) <= FLOATING_POINT_EPSILON) flow_ = 0.0;
    if (abs(flow_ - capacity_) <= FLOATING_POINT_EPSILON) flow_ = capacity_;
    if (!(flow_ >= 0.0)) throw "Flow is negative";
    if (!(flow_ <= capacity_)) throw "Flow exceeds capacity";
}
string Flow_Edge::toString() const {
    stringstream ss;
    ss << v_ << "->" << w_ << " " << flow_ << "/" << capacity_;
    return ss.str();
}
ostream &operator<<(ostream &out, const Flow_Edge &e) {
    out << e.toString();
    return out;
}

bool operator==(const Flow_Edge &e1, const Flow_Edge &e2) {
    // compare edges by flow and capacity
    return e1.flow() == e2.flow() && e1.capacity() == e2.capacity();
}

bool operator!=(const Flow_Edge &e1, const Flow_Edge &e2) {
    return !(e1 == e2);
}

bool operator<(const Flow_Edge &e1, const Flow_Edge &e2) {
    // compare edges by flow and capacity
    return e1.flow() < e2.flow() || (e1.flow() == e2.flow() && e1.capacity() < e2.capacity());
}

bool operator<=(const Flow_Edge &e1, const Flow_Edge &e2) {
    return e1.flow() < e2.flow() || (e1.flow() == e2.flow() && e1.capacity() <= e2.capacity());
}

bool operator>(const Flow_Edge &e1, const Flow_Edge &e2) {
    return e2 < e1;
}

bool operator>=(const Flow_Edge &e1, const Flow_Edge &e2) {
    return e2 <= e1;
}

Flow_Edge::Flow_Edge(const Flow_Edge &e) {
    v_ = e.v_;
    w_ = e.w_;
    capacity_ = e.capacity_;
    flow_ = e.flow_;
}

Flow_Edge &Flow_Edge::operator=(const Flow_Edge &e) {
    if (this != &e) {
        v_ = e.v_;
        w_ = e.w_;
        capacity_ = e.capacity_;
        flow_ = e.flow_;
    }
    return *this;
}

Flow_Edge &Flow_Edge::operator=(Flow_Edge &&e) {
    if (this != &e) {
        v_ = e.v_;
        w_ = e.w_;
        capacity_ = e.capacity_;
        flow_ = e.flow_;
    }
    return *this;
}

#endif //GRAPH_CPP_FLOW_EDGE_H
