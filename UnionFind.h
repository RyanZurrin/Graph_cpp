//
// Created by Ryan.Zurrin001 on 1/11/2022.
//

#ifndef GRAPH_CPP_UNIONFIND_H
#define GRAPH_CPP_UNIONFIND_H
#include <bits/stdc++.h>
using namespace std;

template <typename T>
class UnionFind {
    // the number of elemets in the union find
    int count;
    // used to track the sizes of each component
    int *sz;
    // id[i] is the parent of i
    T *id;
    // used to track the number of components
    int numComponents;
public:
    explicit UnionFind(int n);
    UnionFind(vector<vector<T>> &graph, int n, T union_value = 1);
    ~UnionFind();
    int find(T p);
    void unionElements(T p, T q);
    bool connected(T p, T q);
    int getComponentSize(T p);
    int getNumComponents();
    int getCount();
    int getLargestComponentSize();
    T getParent(T p);
    void unionNeighbors(vector<vector<int>> grid,
                        int i,
                        int j,
                        T union_value = 1);
    void print();

};

template<typename T>
UnionFind<T>::UnionFind(int n) {
    if (n < 0) {
        throw std::invalid_argument("UnionFind: n must be >= 0");
    }
    count = n;
    numComponents = n;
    sz = new int[n];
    id = new T[n];
    for (int i = 0; i < n; i++) {
        sz[i] = 1;
        id[i] = i;
    }
}

template<typename T>
UnionFind<T>::UnionFind(vector<vector<T>> &graph, int n, T union_value) {
    // build the union find data structure from the graph
    if (n < 0) {
        throw std::invalid_argument("UnionFind: n must be >= 0");
    }
    count = n;
    numComponents = n;
    sz = new int[n];
    id = new T[n];
    for (int i = 0; i < n; i++) {
        sz[i] = 1;
        id[i] = i;
    }
    unionNeighbors(graph, 0, 0, union_value);
}
template<typename T>
UnionFind<T>::~UnionFind() {
    delete[] sz;
    delete[] id;
}
template<typename T>
int UnionFind<T>::find(T p) {
    auto root = p;
    while (root != id[root]) {
        root = id[root];
    }
    // path compression
    while (p != root) {
        auto temp = id[p];
        id[p] = root;
        p = temp;
    }
    return root;
}
template<typename T>
void UnionFind<T>::unionElements(T p, T q) {
    T pRoot = find(p);
    T qRoot = find(q);
    if (pRoot == qRoot) {
        return;
    }
    // make smaller root point to larger one
    if (sz[pRoot] < sz[qRoot]) {
        id[pRoot] = qRoot;
        sz[qRoot] += sz[pRoot];
    } else { // make larger root point to smaller one
        id[qRoot] = pRoot;
        sz[pRoot] += sz[qRoot];
    }
    numComponents--;
}
template<typename T>
bool UnionFind<T>::connected(T p, T q) {
    return find(p) == find(q);
}
template<typename T>
int UnionFind<T>::getComponentSize(T p) {
    return sz[find(p)];
}
template<typename T>
int UnionFind<T>::getNumComponents() {
    return numComponents;
}
template<typename T>
int UnionFind<T>::getCount() {
    return count;
}
template<typename T>
T UnionFind<T>::getParent(T p) {
    return id[find(p)];
}
template<typename T>
int UnionFind<T>::getLargestComponentSize() {
    int max = 0;
    for (int i = 0; i < count; i++) {
        if (sz[i] > max) {
            max = sz[i];
        }
    }
    return max;
}
template<typename T>
void UnionFind<T>::print() {
    for (int i = 0; i < count; i++) {
        cout << i << ": " << id[i] << " " << sz[i] << endl;
    }
}
/**
 * @brief iterates over the graph and unions the neighbors of each node that is
 * equal to the union_value
 * @tparam T the type of the graph
 * @param grid the graph
 * @param visited  a boolean array that keeps track of whether a node has been
 * @param i row
 * @param j column
 * @param union_value  the value that is used to determine if a node is a neighbor
 */
template<typename T>
void UnionFind<T>::unionNeighbors(vector<vector<int>> grid,
                                  int i ,
                                  int j, T union_value) {
    auto n = grid.size();
    auto m = grid[0].size();
    int dx[] = {-1, 0, 1, 0};
    int dy[] = {0, 1, 0, -1};
    for(int r = 0; r < n; r++) {
        for(int c = 0; c < m; c++) {
            if (grid[r][c] == union_value) {
                for (int d = 0; d < 4; d++) {
                    int x = r + dx[d];
                    int y = c + dy[d];
                    if (x >= 0 && x < n && y >= 0 && y < m) {
                        if (grid[x][y] == union_value) {
                            unionElements(r * m + c, x * m + y);
                        }
                    }
                }
            }
        }
    }
}
int largest_island(vector<vector<int> > grid, int islValue) {
    // return the size of the largest island using a union find data structure to
    // store the connected components
    int n = grid.size();
    int m = grid[0].size();
    int size = n * m;

    // check grid and see if it is all 0, and if so, return 0
    bool all_zeros = true;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            if (grid[i][j] == 1) {
                all_zeros = false;
                break;
            }
        }
    }
    if (all_zeros) {
        return 0;
    }

    UnionFind<int> uf(grid, size, islValue);
    uf.print();
    // print the number of connected components
    cout << "Number of connected components: " << uf.getNumComponents() << endl;
    // print the size of the largest connected component
    cout << "Largest connected component: " << uf.getLargestComponentSize() << endl;
    // print the number of nodes in the graph
    cout << "Number of nodes: " << uf.getCount() << endl;
    return uf.getLargestComponentSize();
}

#endif //GRAPH_CPP_UNIONFIND_H
