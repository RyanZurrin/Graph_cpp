//#include "Graph.h"
#include <iostream>
#include <chrono>
#include "RunTImer.h"
using namespace std;
#include  "DynamicProblems.h"
int main() {
    RunTimer timer(SECONDS);


    // fll the vecor with random numbers between 1 and 10
//
//    vector<vector<int>> boxes = {
//            {2,  1,  2},
//            {3,  2,  3},
//            {2,  2,  8},
//            {2,  3,  4},
//            {2,  2,  1},
//            {4,  4,  5},
//            {5,  3,  6},
//            {6,  2,  7},
//            {1,  2,  4},
//            {4,  5,  3},
//            {6,  4,  8}
//    };
//    for (int i = 0; i < n; i++) {
//        cin >> stones[i];
//    }
//  initiate a vector with 100 elements all with value 0
    vector<int> stoneHeights = {10, 30, 40, 50, 20};
    timer.start();
    long long bsts = minimumCost(stoneHeights, 3);
    timer.stop();
    cout << "cost of frog to jump is  ";
    cout << bsts << endl;
    timer.display();


    timer.display();

    return 0;
}