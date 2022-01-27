//
// Created by Ryan.Zurrin001 on 1/12/2022. updated
//

#ifndef GRAPH_CPP_DYNAMICPROBLEMS_H
#define GRAPH_CPP_DYNAMICPROBLEMS_H
#include <iostream>
#include <vector>
#include <bits/stdc++.h>
using namespace std;
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
// recursive approach
int longestCommonSubsequence(string s1, string s2, int i, int j) {
    if (i == s1.size() or j == s2.size()) {
        return 0;
    }
    if (s1[i] == s2[i]) {
        return 1 + longestCommonSubsequence(s1, s2, i + 1, j + 1);
    }
    int op1 = longestCommonSubsequence(s1, s2, i + 1, j);
    int op2 = longestCommonSubsequence(s1, s2, i, j + 1);
    return max(op1, op2);
}
// top down approach
int longestCommonSubsequenceTD(string s1, string s2, int i, int j, vector<vector<int>>& dp) {
    if (i == s1.size() or j == s2.size()) {
        return 0;
    }
    if (dp[i][j] != -1) {
        return dp[i][j];
    }
    if (s1[i] == s2[i]) {
        return dp[i][j] = 1 + longestCommonSubsequenceTD(s1, s2, i + 1, j + 1, dp);
    }
    int op1 = longestCommonSubsequenceTD(s1, s2, i + 1, j, dp);
    int op2 = longestCommonSubsequenceTD(s1, s2, i, j + 1, dp);
    return dp[i][j] = max(op1, op2);
}

int longestCommonSubsequenceBU(string s1, string s2) {
    int n = s1.size();
    int m = s2.size();
    vector<vector<int>> dp(n + 1, vector<int>(m + 1, 0));
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            if (s1[i - 1] == s2[j - 1]) {
                dp[i][j] = 1 + dp[i-1][j-1];
            } else {
                int op1 = dp[i - 1][j];
                int op2 = dp[i][j - 1];
                dp[i][j] = max(op1, op2);
            }
        }
    }
    // printing the subsequence
    int i = n, j = m;
    vector<char> ans;
    while (i != 0 and j != 0) {
        if (dp[i][j] == dp[i][j - 1]) {
            j--;
        } else if (dp[i][j] == dp[i - 1][j]) {
            i--;
        } else {
            ans.push_back(s1[i-1]);
            i--;
            j--;
        }
    }
    reverse(ans.begin(), ans.end());
    for (auto c : ans) {
        cout << c << " ";
    }

    return dp[n][m];
}
// top down dp approach
int sellingWines(vector<vector<int>>& dp, vector<int>&prices, int L, int R, int y) {
    if (L > R) {
        return 0;
    }
    if (dp[L][R] != 0) {
        return dp[L][R];
    }

    int pl = y*prices[L] + sellingWines(dp, prices, L + 1, R, y + 1);
    int pr = y*prices[R] + sellingWines(dp, prices, L, R - 1, y + 1);
    return dp[L][R] = max(pl, pr);
}

// bottom up dp approach
int sellingWines(vector<int> prices, int n) {
    vector<vector<int>> dp(n, vector<int>(n + 1, 0));
    for (int i = n - 1; i >= 0; i--) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                dp[i][j] = n*prices[i];
            } else if (i <= j) {
                int year = n - (j - i);
                int pl = prices[i] * year + dp[i + 1][j];
                int pr = prices[j] * year + dp[i][j - 1];
                dp[i][j] = max(pl, pr);
            }
        }
    }
    // printing the subsequence
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << dp[i][j] << " ";
        }
        cout << endl;
    }
    return dp[0][n - 1];
}
// recursive call
int countingSubsequences(string a, string b, int m, int n) {
    if ((m==-1) and (n==-1) or n==-1) {
        return 1;
    }
    if (m==-1) {
        return 0;
    }
    if (a[m] == b[n]) {
        return countingSubsequences(a, b, m - 1, n - 1) + countingSubsequences(a, b, m - 1, n);
    } else {
        return countingSubsequences(a, b, m - 1, n);
    }
}
// bottom up approach
int countingSubsequences(string a, string b) {
    int m = a.size();
    int n = b.size();
    vector<vector<int>> dp(m + 1, vector<int>(n + 1, 0));
    for(int i = 0; i<=m; i++){
        dp[i][0] = 1;
    }
    for(int i = 1; i<=m; i++) {
        for (int j = 1; j <= n; j++) {
            if(a[i-1] == b[j-1]) {
                dp[i][j] = dp[i-1][j-1] + dp[i-1][j];
            } else {
                dp[i][j] = dp[i-1][j];
            }
            cout << dp[i][j]<<" ";
        }
        cout << endl;
    }
    return dp[m][n];
}

// recurive approach
int knapsack(vector<int>&weights, vector<int>&values, int capacity, int n) {
    if (n == 0 or capacity == 0) {
        return 0;
    }
    int inc = 0;
    int exc = 0;
    if (weights[n - 1] <= capacity) {
        inc = values[n - 1] + knapsack(weights, values, capacity - weights[n - 1], n - 1);
    }
    exc = knapsack(weights, values, capacity, n - 1);
    return max(inc, exc);
}
// bottom up approach
int knapsackBU(vector<int> wt, vector<int> price, int N, int W) {
    vector<vector<int>> dp(N + 1, vector<int>(W + 1, 0));
    for (int n = 1; n <= N; n++) {
        for(int w = 1; w <= W; w++) {
            int inc = 0, exc = 0;
            if(wt[n-1] <= w) {
                inc = price[n-1] + dp[n-1][w-wt[n-1]];
            }
            exc = dp[n - 1][w];
            dp[n][w] = max(inc, exc);
        }
    }
    return dp[N][W];
}

// Coin Change 2, given a value N and an integer vector COINS representing coins of
// different denominations. COnsider you have infinite supply of each coin, your
// task is to find the total number of combinations of these coins that make a sum of N.
// If that amount of money cannot be made up by any combination of coins, return  0;
long long int coinChange(vector<int>& coins, int n) {
    vector<vector<int>> dp(n + 1, vector<int>(coins.size(), 0));
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j < coins.size(); j++) {
            if (i == 0) {
                dp[i][j] = 1;
            } else if (i < coins[j]) {
                dp[i][j] = dp[i][j - 1];
            } else {
                dp[i][j] = dp[i][j - 1] + dp[i - coins[j]][j];
            }
        }
    }
    return dp[n][coins.size() - 1];
}

// game of wits helper function
int gameOfWitsUtil(string s, int i, int j, vector<vector<int>>& dp) {
    if(i+1==j){
        if(s[i]=='O' && s[j]=='O'){
            return 2;
        }
        else if(s[i]=='H'&& s[j]=='H'){
            return -3;
        }
        else{
            return -1;
        }
    }
    if(i==j){
        if(s[i]=='O'){
            return 1;
        }
        else{
            return -2;
        }
    }
    if(dp[i][j] != 0) return dp[i][j];


    if(s[i]=='H' && s[j]=='H'){
        return dp[i][j] = -1* (j-i+1) - 1;
    }
    else if(s[i]=='H' && s[j]=='O'){
        int op1=std::numeric_limits<int>::max();
        if(s[j-1]=='H'){
            op1 = gameOfWitsUtil(s,i,j-2,dp);
        }
        int op2 = gameOfWitsUtil(s,i+1,j-1,dp);
        return dp[i][j] =  min(op1,op2);
    }
    else if(s[i]=='O' && s[j]=='H'){
        int op1=std::numeric_limits<int>::max();
        if(s[i+1]=='H'){
            op1 = gameOfWitsUtil(s,i+2,j,dp);
        }
        int op2 = gameOfWitsUtil(s,i+1,j-1,dp);
        return dp[i][j] =  min(op1,op2);
    }
    else{

        if(s[i+1]=='O'){
            return j-i+1;
        }
        else if(s[j-1]=='O'){

            return j-i+1;
        }
        else {
            int op1=gameOfWitsUtil(s,i+2,j,dp);
            int op2=gameOfWitsUtil(s,i,j-2,dp);
            return dp[i][j] =  max(op1,op2);
        }
    }
}

// Oswald and Henry are playing a game with alternating turns. Oswald always plays first.
// The game starts with all game cards arranged in a single row. The cards are of two
// types 'O' and 'H'. During Oswald's turn, he must choose and remove an 'O' card that
// is either the leftmost or rightmost card remaining. During Henry's turn, he must choose
// and remove an 'H' card that is either the leftmost or rightmost card remaining. IF
// at any point one of the players does not have a legal move (possibly because he has
// no cards left), he loses the game and the other player is awarded 1 point plus 1
// addidional point for  each card that remians on the board.
pair<char,int> GameOfWits(string s){
    int n = s.size();
    vector<vector<int>> dp(n,vector<int>(n,0));
    int res = gameOfWitsUtil(s,0,n-1,dp);
    if (res > 0) {
        return {'O',abs(res)};
    } else {
        return {'H',abs(res)};
    }
}
long long smokeMixUtil(vector<int> vec, int i, int j, vector<vector<long long>>&dp) {
    if (i == j) {
        return 0.0;
    }
    if (dp[i][j] != -1) {
        return dp[i][j];
    }
    long long res = std::numeric_limits<long long>::max();
    for (int k = i; k < j; k++) {
        auto result = smokeMixUtil(vec, i, k, dp) + smokeMixUtil(vec, k + 1, j, dp);
        result += ( ( (vec[k] - (i > 0 ? vec[i - 1] : 0) ) % 100 ) * ((vec[j] - vec[k]) % 100) );
        res = min(res, result);
    }
    return dp[i][j] = res;
}
// Mixtures-SPOJ
// Harry Potter has n mixturs in front of him, arrainged in a row. Each mixture has one
// of 100 different colors (colors have numbers from 0 to 99). He wants to mix all
// these mixtures together. At each step, he is going to take two mixtures that stand
// next to each other and mix them together, and put the resulting mixture in thier place.
// When mixing two mixtures of colors a and b, the resulting misture will ahve the color
// (a+b)%100. Also, there will be some smoke in the process. THe amount of smoke is
// generated when mixing two mixtures of colors a and b and is a*b. FInd out
// what is the minimum amount of smoke that will be generated if Harry Potter mixes all
// the mixtures together.
long long minimumSmoke(vector<int> v){
    int n = v.size();
    vector<vector<long long>> dp(n,vector<long long>(n,-1));
    for (int i = 1; i < n; i++) {
        v[i] += v[i - 1];
    }
    auto res = smokeMixUtil(v, 0, n - 1, dp);
    return res;
}

// Busy Life
//  you are actually very busy person. You have a long list of activities that you
//  want to do. Your aim is to finish as many activities as possible. In the given
//  arrangment if you go to one activity  and the tiems conflict with the next activity
//  you will have to skip that activity. FInd the maximum number of activities that
//  you can finish.
int countActivites(vector<pair<int,int> > activities){
    sort(activities.begin(),activities.end() ,[](pair<int,int> a,pair<int,int> b){
        return a.second < b.second;
    });
    int count = 1;
    int completed = activities[0].second;
    for (int i = 1; i < activities.size(); i++) {
        if(activities[i].first >= completed){
            count++;
            completed = activities[i].second;
        }
    }
    return count;
}

int badness(vector<pair<string,int> > teams){
    int n = teams.size();
    vector<int> dp(n+1,0);
    for (int i = 0; i < n; i++) {
        dp[teams[i].second]++;
    }
    int res = 1;
    int sum = 0;
    for (int i = 1; i <= n; i++) {
        while (dp[i]) {
            sum += abs(res - i);
            dp[i]--;
            res++;
        }
    }
    return sum;
}

int minDistance(string a, string b, int na, int nb, vector<vector<int>>&dp){
    if (na == 0) {
        return nb;
    }
    if (nb == 0) {
        return na;
    }
    if (dp[na][nb] != -1) {
        return dp[na][nb];
    }
    int res = std::numeric_limits<int>::max();
    if (a[na - 1] == b[nb - 1]) {
        res = minDistance(a, b, na - 1, nb - 1, dp);
    } else {
        int op1 = minDistance(a, b, na - 1, nb, dp);
        int op2 = minDistance(a, b, na, nb - 1, dp);
        int ob3 = minDistance(a, b, na - 1, nb - 1, dp);
        res = min(op1, min(op2, ob3)) + 1;
    }
    return dp[na][nb] = res;
}
// Given two strings str1 and str2, find the minimum number of edits (operations)
// required to convert str1 equal to str2. The following operations are allowed:
// 1. Insert a character
// 2. Remove a character
// 3. Replace a character
int editDistance(string str1, string str2){
    int n = str1.size();
    int m = str2.size();
    vector<vector<int>> dp(n+1,vector<int>(m+1,-1));
    int res = minDistance(str1, str2, n, m, dp);
    return res;
}

int hasMatch(string s, string p,int n, int m, vector<vector<int>>&dp){
    if (n == 0 && m == 0) {
        return 1;
    }
    if (m == 0) {
        return 0;
    }
    if (n == 0) {
        while(m > 0) {
            if (p[m-1] != '*') {
                return 0;
            }
            m--;
        }
        return 1;
    }
    if (dp[n][m] != -1) {
        return dp[n][m];
    }
    int res = 0;
    if (dp[n][m] != -1) {
        return dp[n][m];
    }
    if (s[n-1] == p[m-1] || p[m-1] == '?') {
        res = hasMatch(s, p, n-1, m-1, dp);
    } else if (p[m-1] == '*') {
        res = hasMatch(s, p, n-1, m, dp) || hasMatch(s, p, n, m-1, dp);
    } else {
        res = 0;
    }
    return dp[n][m] = res;
}


bool wildcardPatternMatching(string s, string p){
    int n = s.size();
    int m = p.size();
    vector<vector<int>> dp(n+1,vector<int>(m+1,-1));
    int res = hasMatch(s, p, n, m, dp);
    if (res) {
        return true;
    }
    return false;
}

int palindromicPartitioningUtil(string s){
    // min cuts to make s palindromic
    int n = s.length();
    vector<vector<bool> > isPalin(n+1,vector<bool>(n,false));

    for(int i=0;i<n;i++){
        isPalin[i][i] = true;
    }

    //2d dp palindromic grid for helper
    // tell whether a string i...j is a palindrome or not
    for(int len=2;len<=n;len++){
        for(int i=0;i<=n-len;i++){
            //substring i to j
            int j = i + len - 1;
            if(len==2){
                isPalin[i][j] = (s[i]==s[j]);
            }
            else{
                isPalin[i][j] = (s[i]==s[j] and isPalin[i+1][j-1]);
            }
        }
    }

    //min cut logic
    vector<int> cuts(n+1,INT_MAX);

    for(int i=0;i<n;i++){
        if(isPalin[0][i]){
            cuts[i] = 0;
        }
        else{
            cuts[i]  = cuts[i-1] + 1;
            for(int j=1;j<i;j++){
                if(isPalin[j][i] and cuts[j-1] + 1 < cuts[i]){
                    cuts[i] = cuts[j-1] + 1;
                }
            }
        }
    }
    return cuts[n-1];
}

int palindromicPartitioning(string str) {
    int res = palindromicPartitioningUtil(str);
    return res;
}

#endif //GRAPH_CPP_DYNAMICPROBLEMS_H
