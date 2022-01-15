//
// Created by Ryan.Zurrin001 on 1/12/2022.
//

#ifndef GRAPH_CPP_DYNAMICPROBLEMS_H
#define GRAPH_CPP_DYNAMICPROBLEMS_H
#include <iostream>
#include <vector>
// O(k^n)
long long countWays(int n, int k) {
    if (n == 0) {
        return 1;
    }
    if (n < 0 || k < 0) {
        return 0;
    }
    auto ways = 0;
    for (int jump = 1; jump <= k; jump++) {
        ways += countWays(n - jump, k);
    }
    return ways;
}

// top down approach O(nk)
long long countWaysTD(int n, int k, int *dp) {
    if (n == 0) {
        return 1;
    }
    if (n < 0) {
        return 0;
    }
    if (dp[n] != 0) {
        return dp[n];
    }
    auto ways = 0;
    for (int jump = 1; jump <= k; jump++) {
        ways += countWaysTD(n - jump, k, dp);
    }
    return dp[n] = ways;
}

// bottom up approach O(nk)
long long countWaysBU(int n, int k) {
    vector<long long> dp(n+1, 0);
    dp[0] = 1;
    for (int i = 1; i <= n; i++) {
        for (int jump = 1; jump <= k; jump++) {
            if (i - jump >= 0) {
                dp[i] += dp[i - jump];
            }
        }
    }
    return dp[n];
}

long long countWaysOptBU(int n, int k) {
    vector<long long> dp(n+1, 0);
    dp[0] = dp[1] = 1;
    for (int i = 2; i <= n; i++) {
        for (int jump = 1; jump <= k; jump++) {
            if (i - jump >= 0) {
                dp[i] += dp[i - jump];
            }
        }
    }
    return dp[n];
}



#endif //GRAPH_CPP_DYNAMICPROBLEMS_H
