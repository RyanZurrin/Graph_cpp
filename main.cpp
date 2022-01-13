#include "Graph.h"
#include "MaxFlowGraph.h"

int main() {
    int V = 5;
    int E = 10;
    int s = 0;
    int t = V-1;
    FlowNetwork G(V, E);
    cout << G << endl;

    FordFulkerson maxFlow(G, s, t);
    cout << "Max flow from " << s << " to " << t << " is:\n";
    maxFlow.printMaxFlow();
    // print min cut
    cout << "Min cut:\n";
    maxFlow.printMinCut();

    // print value
    cout << "Value: " << maxFlow.getValue() << endl;


    return 0;
}