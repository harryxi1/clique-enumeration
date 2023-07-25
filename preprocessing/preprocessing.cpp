#include <iostream>
#include <vector>
#include <algorithm>

#include <fstream>
#include <string>
#include <sstream>

#include <chrono>

using namespace std;

#define FILENAME "../data/web_Stanford_edges.csv"
#define n 281903

// lighter weight standalone preprocessing

struct node {
int label;
int deg;
int colour;

std::vector<int> neighbours; // N(v) in undirected graph

};

int getData(vector<node> &vertices){

    string filename = FILENAME;
    string line;

    ifstream data;
    data.open(filename);
    int maxDegree = 0;

    while (getline(data, line)) {
        istringstream iss(line);
        string token;
        getline(iss, token, ',');
        int u = atoi(token.c_str());

        getline(iss, token, ',');
        int v = atoi(token.c_str());
        if (u != v) {
            vertices[u].deg++;
            vertices[u].neighbours.push_back(v);

            vertices[v].deg++;
            vertices[v].neighbours.push_back(u);
            maxDegree = max({maxDegree, vertices[u].deg, vertices[v].deg++});
        }
    }
    return maxDegree;
}

vector<int> degeneracyOrdering(vector<node> G, int maxDegree) {
    int k = 0;
    vector<bool> marked (n, 0); // 1 if in L, otherwise 0
    vector<int> d (n, 0); // dv 

    vector<int> L; L.reserve(n);
    vector<vector<int>> D (maxDegree);

    for (node v : G) {
        D[v.deg].push_back(v.label); 
        d[v.label] = v.deg;
    }

    for (int i = 0; i < n; i++) {
        if (i % 25000 == 0) {cout << i << "\n";}
        int currBucket = 0;
        while (D[currBucket].empty()) {currBucket++;}
        k = max(k, currBucket);
        node v = G[D[currBucket].back()]; // picks last vertex in D[currBucket]
        D[currBucket].pop_back(); // remove it from the bucket
        L.push_back(v.label); marked[v.label] = 1; // add it to the BACK of L 
        for (int w : v.neighbours) {
            if (marked[w] == 0) {
                D[d[w]].erase(remove(D[d[w]].begin(), D[d[w]].end(), w), 
                              D[d[w]].end());
                D[d[w]-1].push_back(w);
                d[w]--;
            }
        }
    }
    return L;
}

int greedyColouring(vector<node> &G, vector<int> ordering) {
    G[ordering[n-1]].colour = 0;
    int r = 1;
    for (int v = G.size()-2; v >= 0; v--) {
        if (v % 25000 == 0) {cout << v << "\n";}
        vector<bool> available(G.size(), 1);
        for (int w : G[ordering[v]].neighbours) {
            if (G[w].colour != -1) {
                available[G[w].colour] = 0;
            }
        }
        for (int c = 0; c < G.size(); c++) {
            if (available[c] == 1) {
                G[ordering[v]].colour = c;
                r = max(r, c);
                break;
            }
        }
    }
    return r;
}

int main() {
    vector<node> G; G.reserve(n);
    for (int v = 0; v < n; v++) {
        node curr;
        curr.label = v; curr.deg = 0; curr.colour = -1;
        G.push_back(curr);
    }; 
    int maxDegree = getData(G); cout << "Data Retrieved\n";
    auto t1 = chrono::high_resolution_clock::now();
    vector<int> degOrdering = degeneracyOrdering(G, maxDegree);
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Found Degeneracy Ordering in (" << 
    chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms)\n";


    greedyColouring(G, degOrdering);
    vector<int> colours;
    for (auto v : G) {
        colours.push_back(v.colour);
    }

    ofstream degOutFile("degStanford.txt");
    for (const auto &e : degOrdering) degOutFile << e << "\n";

    ofstream degOutFile("colourStanford.txt");
    for (const auto &e : colours) degOutFile << e << "\n";
}