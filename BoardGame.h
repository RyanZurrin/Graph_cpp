//
// Created by Ryan.Zurrin001 on 1/10/2022.
//

#ifndef GRAPH_CPP_BOARDGAME_H
#define GRAPH_CPP_BOARDGAME_H
#include <iostream>
#include <unordered_map>
using namespace std;

class Node {
public:
    char s;
    unordered_map<char, Node*> children;
    string word;
    bool isTerminal;
    explicit Node(char s) {
        this->s = s;
        isTerminal = false;
        word = "";
    }
};

class Trie {
public:
    Node* root;
    Trie() {
        root = new Node(' ');
    }
    void addWord(const string& word) const {
        Node* curr = root;
        for (auto ch : word) {
            if (curr->children.count(ch) == 0) {
                curr->children[ch] = new Node(ch);
            }
            curr = curr->children[ch];
        }
        curr->isTerminal = true;
        curr->word = word;
    }

    [[maybe_unused]] bool search(const string& word) const {
        Node* curr = root;
        for (auto ch : word) {
            if (curr->children.count(ch) == 0) {
                return false;
            }
            curr = curr->children[ch];
        }
        return curr->isTerminal;
    }
};


void dfs(vector<vector<char>> board, Node *node, int i, int j,
         vector<vector<bool>> visited, unordered_set<string>& output) {
    // Base case
    char ch = board[i][j];
    if (node->children.count(ch) == 0) {
        return;
    }
    visited[i][j] = true;
    node = node->children[ch];
    if (node->isTerminal) {
        output.insert(node->word);
    }

    // make 8 direction arrays
    int dx[] = {0, 1, 1, 1, 0, -1, -1, -1};
    int dy[] = {-1, -1, 0, 1, 1, 1, 0, -1};
    for (int k = 0; k < 8; k++) {
        int x = i + dx[k];
        int y = j + dy[k];
        if (x >= 0 and x < board.size() and y >= 0 and y < board[0].size() &&
            !visited[x][y]) {
            dfs(board, node, x, y, visited, output);
        }
    }
    visited[i][j] = false;
}

[[maybe_unused]] void playBoggle(vector<string> words, vector<vector<char>> board, int m, int n) {
    Trie trie;
    for (auto word : words) {
        trie.addWord(word);
    }
    vector<vector<bool>> visited(m, vector<bool>(n, false));
    unordered_set<string> foundWords;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            dfs(board, trie.root, i, j, visited, foundWords);
        }
    }
    for (auto word : foundWords) {
        cout << word << endl;
    }
}

/* Driver program to test above code
 * const int M = 19;
    const int N = 5;

    vector<string> words = {"SNAKE", "FOR", "QUE", "SNACK",
                      "SNACKS", "GO", "TUNES", "CAT", "BUZZ", "HIGH",
                      "CODE", "QUIZ", "DIDO", "CODING", "MONKEY", "BUFF",
                      "RICH", "ENIGMA", "GOD", "SMELT", "HACK", "ZOO",
                      "ZOOM", "RAKE", "MELT", "PROGRAMMER"};
    // BUILD A 8 BY 8 GRAPH
    vector<vector<char>> boggle = {
            {'S', 'N', 'A', 'K', 'E'},
            {'F', 'O', 'R', 'Q', 'U'},
            {'E', 'Z', 'X', 'S', 'U'},
            {'P', 'N', 'A', 'C', 'K'},
            {'S', 'R', 'A', 'C', 'K'},
            {'G', 'O', 'A', 'M', 'M'},
            {'E', 'R', 'E', 'E', 'B'},
            {'C', 'O', 'D', 'R', 'U'},
            {'Q', 'U', 'I', 'Z', 'Z'},
            {'D', 'D', 'D', 'L', 'O'},
            {'C', 'O', 'I', 'N', 'G'},
            {'M', 'O', 'N', 'M', 'Y'},
            {'B', 'U', 'G', 'F', 'A'},
            {'R', 'I', 'C', 'H', 'E'},
            {'E', 'N', 'I', 'D', 'M'},
            {'A', 'G', 'O', 'Y', 'E'},
            {'S', 'M', 'E', 'L', 'T'},
            {'H', 'A', 'C', 'K', 'N'},
            {'O', 'Z', 'O', 'O', 'M'}
    };
    playBoggle(words, boggle, M, N);
 */




#endif //GRAPH_CPP_BOARDGAME_H
