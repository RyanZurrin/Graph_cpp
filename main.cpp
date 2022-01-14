//#include "Graph.h"

#include "Ford_Fulkerson.h"
int main() {

//    Flow_Network G(1000, 10000);
//    cout << G << endl;
    int V, E, s, t;
    V = 10;
    E = 20;
    s = 0;
    t = V-1;

    Flow_Network f(V, E);
    cout << f << endl;

    Ford_Fulkerson ff(f, s, t);
    ff.ford_fulkerson_algorithm();




    return 0;
}