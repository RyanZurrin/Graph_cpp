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
// bottom up approach O(nk)
int minNumberOfCoinsForChange(int n, vector<int> &c) {
    vector<int> dp(n+1, 0);
    dp[0] = 0;
    for (int i = 1; i <= n; i++) {
        dp[i] = std::numeric_limits<int>::max();
        for(int coin : c) {
            if (i - coin >= 0 && dp[i - coin] != std::numeric_limits<int>::max()) {
                dp[i] = min(dp[i], dp[i - coin] + 1);
            }
        }
    }
    return dp[n] == std::numeric_limits<int>::max() ? -1 : dp[n];
}

// recursive approach O(nk)
int maxProfit(int prices[], int n) {
    if (n <= 0) {
        return 0;
    }
    int maxProfit_ = std::numeric_limits<int>::min();
    for(int i = 0; i < n; i++) {
        int cut = i + 1;
        int current_ans = prices[i] + maxProfit(prices, n - cut);
        maxProfit_ = max(maxProfit_, current_ans);
    }
    return maxProfit_;
}
// bottom up approach O(nk)
int maxProfitBU(int prices[], int n) {
    vector<int> dp(n+1, 0);
    dp[0] = 0;
    for (int len = 1; len <= n; len++) {
        int ans = std::numeric_limits<int>::min();
        for (int i = 0; i < len; i++) {
            int cut = i + 1;
            int current_ans = prices[i] + dp[len - cut];
            ans = max(ans, current_ans);
        }
        dp[len] = ans;
    }
    return dp[n];
}

int min_jumpsBU(vector<int> &arr) {
    if (arr.size() == 0) {
        return 0;
    }
    vector<int> dp(arr.size(), 0);
    dp[0] = 0;
    for (int i = 1; i < arr.size(); i++) {
        dp[i] = std::numeric_limits<int>::max();
        for (int j = 0; j < i; j++) {
            if (arr[j] + j >= i) {
                dp[i] = min(dp[i], dp[j] + 1);
            }
        }
    }
    return dp[arr.size() - 1];
}

int min_jumpsTD(vector<int> arr, int n, vector<int>& dp, int i =0) {
    if (i == n - 1) {
        return 0;
    }
    if (dp[i] != 0) {
        return dp[i];
    }
    if (i >= n) {
        return std::numeric_limits<int>::max();
    }
    int steps = std::numeric_limits<int>::max();
    int maxJump = arr[i];
    for (int jump = 1; jump <= maxJump; jump++) {
        int next = i + jump;
        int sub_prob = min_jumpsTD(arr, n, dp, next);
        if (sub_prob != std::numeric_limits<int>::max()) {
            steps = min(steps, sub_prob + 1);
        }
    }
    return dp[i] = steps;
}



int getMaxNonAdjacentSum(vector<int> &arr) {
    int n = arr.size();
    if (n == 1) {
        return max(arr[0], 0);
    } else if (n == 2) {
        return max(0 , max(arr[0], arr[1]));
    }
    vector<int> dp(n+1, 0);
    dp[0] = max(0, arr[0]);
    dp[1] = max(0, max(arr[0], arr[1]));
    for (int i = 2; i < n; i++) {
        dp[i] = max(dp[i - 1], dp[i - 2] + arr[i]);
    }
    return dp[n - 1];
}

int longestIncreasingSubsequence(vector<int> &arr) {
    int n = arr.size();
    vector<int> dp(n, 1);
    int len = 1;
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            //cout << "i: " << i << " j: " << j << endl;
            if (arr[i] > arr[j]) {
                //cout << "a[i] > a[j] " << endl;
                //cout << "dp[j] + 1: " << dp[j]+1 << ", dp[i]: " << dp[i] << endl;
                dp[i] = max(dp[i], dp[j] + 1);
                len = max(len, dp[i]);
                //cout << "dp[i]: " << dp[i] << ", len: " << len << endl;
            }
        }
    }
    return len;
}
bool compareBoxes(vector<int> &a, vector<int> &b) {
    return a[2] > b[2];
}
bool canPlace(vector<int> &a, vector<int>&b) {
    return a[0] > b[0] && a[1] > b[1] && a[2] > b[2];
}

int boxStacking(vector<vector<int>> boxes) {
    sort(boxes.begin(), boxes.end(), compareBoxes);
    int n = boxes.size();
    vector<int> dp(n+1, 0);

    for (int i = 0; i < n; i++) {
        dp[i] = boxes[i][2];
    }
    for (int i = 1; i < n; i++) {
        for (int j = 0; j < i; j++) {
            if (canPlace(boxes[j], boxes[i])) {
                int cur_height = boxes[i][2];
                if (dp[j] + cur_height > dp[i]) {
                    dp[i] = dp[j] + cur_height;
                }
            }
        }
    }
    int max_height = 0;
    for (int i = 0; i < n; i++) {
        cout << dp[i] << " ";
        max_height = max(max_height, dp[i]);
    }
    return max_height;
}
//______________________________________________________________________________
// this is not possilbe to calculate n = 35, at around n = 22 it startes taking
// in the multiple miniute range.
unsigned long long int countBST(int n) {
    if (n == 0 || n == 1) {
        return 1;
    }
    unsigned long long int ans = 0;
    for (int i = 1; i <= n; i++) {
        unsigned long long int x = countBST(i - 1);
        unsigned long long int y = countBST(n - i);
        ans += x * y;
    }
    return ans;
} // end recursive apporach
// this takes less than 1 second to calculate n = 35
unsigned long long int countBST_TD(int n, vector<unsigned long long int>& dp) {
    if (n == 0 || n == 1) {
        return 1;
    }
    if (dp[n] != 0) {
        return dp[n];
    }
    unsigned long long int ans = 0;
    for (int i = 1; i <= n; i++) {
        unsigned long long int x = countBST_TD(i - 1, dp);
        unsigned long long int y = countBST_TD(n - i, dp);
        ans += x * y;
    }
    return dp[n] = ans;
}
// this is the fastes of the three by far
unsigned long long int countBST_BU(int n) {
    vector<unsigned long long int> dp(n+1, 0);
    dp[0] = dp[1] = 1;
    for (int i = 2; i <= n; i++) {
        for (int j = 1; j <= i; j++) {
            dp[i] += dp[j - 1] * dp[i - j];
        }
    }
    return dp[n];
}
//______________________________________________________________________________

int getMinCostBU(vector<int> stones) {
    int n = stones.size();
    vector<int> dp(n, 0);
    dp[1] = abs(stones[1] - stones[0]);
    for (int i = 2; i < n; i++ ) {
        int op1 = abs(stones[i] - stones[i - 1]) + dp[i - 1];
        int op2 = abs(stones[i] - stones[i - 2]) + dp[i - 2];
        dp[i] = min(op1, op2);
    }
    return dp[n - 1];
}

// there are N stones, for each i (1 <= i <= N), the height of stone i is h[i]. There
// is a frog who is initially on stone 1. He will repeat the following action some number
// of times to reach stone N: if the frog is currently on stone i, jump to one of the
// following: stone i + 1, i + 2,..., i+k. Here a cost of [hi-hj] is icurred, where
// j is the stone to land on.
// find the minimum possible total cost incurred before the from reaches stone N.
long long minimumCost(vector<int> stones, int k){
    //Write your code here. Do not modify the function or the parameters. Use a helper function if needed.
    int n = stones.size();
    vector<long long> dp(n);
    for (int i = 1; i < n; i++) {
        dp[i] = std::numeric_limits<long long>::max();
    }
    // print the dp array
//    for (int i = 0; i < n; i++) {
//        cout << dp[i] << " ";
//    }
//    cout << endl;
//    // print out the stones array
//    for (int i = 0; i < n; i++) {
//        cout << stones[i] << " ";
//    }
//    cout << endl;

    for (int i = 0; i < n; i++) {
        for (int j = i+1; j <= i+k; j++) {
            if (j < n) {
                //cout <<"min: " << dp[i] << " + " << abs(stones[i] - stones[j]) << " = ";
                dp[j] = min(dp[j], dp[i] + abs(stones[j] - stones[i]));
                //cout << dp[j] << " ";
            }
        }
    }
    return dp[n - 1];
}

#endif //GRAPH_CPP_DYNAMICPROBLEMS_H
