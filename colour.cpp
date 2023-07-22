#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
#include <random>
#include <chrono>

#include "utils.h"
#include "graphutils.h"
#include "pivoter.h"
#include "suffix.h"

using namespace std;

random_device rd;
mt19937 gen(rd());

#define FILENAME "data/web_Stanford_edges.csv"
#define n 281903
#define m 2312497

void getData(vector<node> &vertices){

    string filename = FILENAME;
    string line;

    ifstream data;
    data.open(filename);

    while (getline(data, line)) {
        istringstream iss(line);
        string token;
        int counter = 0;
        getline(iss, token, ',');
        int u = atoi(token.c_str());

        getline(iss, token, ',');
        int v = atoi(token.c_str());
        counter++;
        if (counter % 10000 == 0) {
            cout << counter << endl;
        }
        if (u != v) {
            vertices[u].deg++;
            vertices[u].neighbours.push_back(v);

            vertices[v].deg++;
            vertices[v].neighbours.push_back(u);
        }
    }
}

int greedyColouring(vector<node> &vertices, vector<int> degOrdering) {
    vector<int> ordering (degOrdering.rbegin(), degOrdering.rend());
    vector<int> available(vertices.size(), 1);
    vertices[ordering[0]].colour = 0;
    int r = 1;
    for (int v = 1; v < n; v++) {
        for (int i = 0; i < vertices[ordering[v]].neighbours.size(); i++) {
            node curr_neigh = vertices[vertices[ordering[v]].neighbours[i]];
            if (curr_neigh.colour < n + 1 && available[curr_neigh.colour] == 1) {
                available[curr_neigh.colour] = 0;
            }
        }
        int i = 0;
        for (i = 0; i < n; i++)
            if (available[i] == 1)
                break;
            else {available[i] = 1;}

        vertices[ordering[v]].colour = i;
        if (i > r) {r = i;}
        for (int j = i+1; j < n; j++) {
            if (available[j] == 0) {available[j] = 1;}
        }
    }
    return r;
}

void setDAGNeighbourhoods(vector<node> &vertices, vector<int> ordering) {
    unordered_map<int, int> ordering_index_map = {};
    for (int i = 0; i < n; i++) {
        ordering_index_map[ordering[i]] = i;
    }
    for (int v = 0; v < n; v++) {
        for (int i = 0; i < vertices[v].neighbours.size(); i++) {
            if (ordering_index_map[v] < ordering_index_map[vertices[v].neighbours[i]]) {
                vertices[v].outNeighbours.push_back(vertices[v].neighbours[i]);
            }
        }
        // sort to save on intersection/difference costs
        sort(vertices[v].neighbours.begin(), vertices[v].neighbours.end());
        sort(vertices[v].outNeighbours.begin(), vertices[v].outNeighbours.end());

        // remove potential duplicates
        //vertices[v].neighbours.erase(unique(vertices[v].neighbours.begin(), vertices[v].neighbours.end() ), vertices[v].neighbours.end());
        //vertices[v].outNeighbours.erase(unique(vertices[v].outNeighbours.begin(), vertices[v].outNeighbours.end() ), vertices[v].outNeighbours.end());
    }
}

void setColourDAGNeighbourhoods(vector<node> &vertices, vector<int> colourOrdering) {
    unordered_map<int, int> ordering_index_map = {};
    for (int i = 0; i < n; i++) {
        ordering_index_map[colourOrdering[i]] = i;
    }

    for (int v = 0; v < n; v++) {
        for (int i = 0; i < vertices[v].neighbours.size(); i++) {
            if (ordering_index_map[v] < ordering_index_map[vertices[v].neighbours[i]]) {
                vertices[v].colourOutneighbours.push_back(vertices[v].neighbours[i]);
            }
        }
        sort(vertices[v].colourOutneighbours.begin(), vertices[v].colourOutneighbours.end());
    }
}

// bool partitioner0(int k, node v, vector<node> vertices) {
//     vector<int> neighbours = v.outneighbours;
//     int neighbourhood_size = neighbours.size();
//     float n_neighbour = static_cast<float>(neighbourhood_size);
//     if (n_neighbour == 0) {
//         return 0;
//     }
//     float degsum = 0;

//     for (int u = 0; u < n_neighbour; u++) {
//         int neighbour_degree = vertices[neighbours[u]].neighbours.size();
//         vector<int> u_neighbours = vertices[neighbours[u]].neighbours;
//         int i = 0, j = 0;
//         while (j < neighbourhood_size && i < neighbour_degree) {
//             if (u_neighbours[i] < neighbours[j]) {i++;}
//             else if (neighbours[j] < u_neighbours[i]) {j++;}
//             else {degsum++; i++; j++;}
//         }
//     }
//     return degsum/n_neighbour > k;
// }

bool partitioner(int k, node v, vector<node> outNeighbours) {
    float neighbourhoodSize = outNeighbours.size();
    if (neighbourhoodSize == 0) {return 0;}
    float degsum = 0;
    vector<node> inducedSubgraph = induceSubgraph(outNeighbours);
    for (node v : inducedSubgraph) {degsum = degsum + v.neighbours.size();}
    return degsum/neighbourhoodSize > k;
}

vector<vector<int>> partitionVertices(int k, vector<node> vertices) {
    // returns V partititon into "non dense" and "dense" (S) by outneighbourhood
    vector<vector<int>> result = {{}, {}};
    for (int v = 0; v < n; v++) {
        vector <node> outNeighbours;
        vector <int> outNeighboursIndex = vertices[v].outNeighbours;
        for (int u = 0; u < outNeighboursIndex.size(); u++) {
            outNeighbours.push_back(vertices[outNeighboursIndex[u]]);
        }
        if (partitioner(k, vertices[v], outNeighbours)) {result[1].push_back(v);}
        else {result[0].push_back(v);}
    }
    return result;
}

auto colourCMP(unordered_map<int, int> colourMap) {
    return [&colourMap](int lhs, int rhs) {
        if (colourMap[lhs] == colourMap[rhs]) {
            return lhs < rhs;
        }
        else {return colourMap[lhs] < colourMap[rhs];}
    };
}

vector<int> getColourOrder(vector<node> vertices) {
    vector<int> colourOrder;
    unordered_map<int, int> colourMap = {};
    for (int i = 0; i < n; i++) {
        colourOrder.push_back(i);
        colourMap[i] = vertices[i].colour;
    }
    sort(colourOrder.begin(), colourOrder.end(), colourCMP(colourMap));
    return colourOrder;
}

vector<vector<int>> DPPathCount(vector<node> vertices, unordered_map<int, int> labelIndexMap, int k) {
    vector<vector<int>> H(vertices.size(), vector<int>(k, 0));
    for (int i = 0; i < vertices.size(); i++) {H[i][0] = 1;}

    for (int j = 1; j < k; j++) {
        for (int i = 0; i < vertices.size(); i++) {
            for (int x = 0; x < vertices[i].colourOutneighbours.size(); x++) {
                int currNeighbourIndex = labelIndexMap[vertices[i].colourOutneighbours[x]];
                node currNeighbour = vertices[currNeighbourIndex];
                H[i][j] = H[i][j] + H[labelIndexMap[currNeighbour.label]][j-1];
            }
        }
    }
    return H;
}

vector<int> DPPathSampler(vector<node> vertices, unordered_map<int, int> labelIndexMap, int k) {
    vector<vector<int>> H = DPPathCount(vertices, labelIndexMap, k);
    vector<int> R = {};
    vector<int> Q = {};
    for (int i = 0; i < vertices.size(); i++) {Q.push_back(vertices[i].label);}
    for (int i = 0; i < k; i++) {
        int cnt = 0;
        for (int u = 0; u < Q.size(); u++) {cnt = cnt + H[labelIndexMap[Q[u]]][k - i - 1];}

        if (cnt == 0) {
            return {};
        }

        discrete_distribution<> distr(Q.begin(), Q.end());
        vector<double> p = distr.probabilities();
        for (int u = 0; u < Q.size(); u++) {
            p[u] = H[labelIndexMap[Q[u]]][k - i - 1]/cnt;
        }
        int u = distr(gen);
        R.push_back(Q[u]); 
        Q = vertices[labelIndexMap[Q[u]]].colourOutneighbours;
    }
    return R;
}

float estimateClique(vector<node> G, vector<int> S, int k, int t, int r) {
    unordered_map<int, vector<vector<int>>> F = {};
    int cntKCol = 0;
    for (auto v : S) {
        // build induced subgraph w/o function
        vector<node> vertexSubset; 
        vertexSubset.reserve(G[v].outNeighbours.size());
        unordered_map<int, int> labelIndexMap = {};

        for (int i = 0; i < G[v].outNeighbours.size(); i++) {
            node vertexCopy = G[G[v].outNeighbours[i]];
            labelIndexMap[vertexCopy.label] = i;
            vertexSubset.push_back(vertexCopy);
        }
        F[v] = DPPathCount(induceSubgraph(vertexSubset), labelIndexMap, k-1);
    }

    unordered_map<int, int> paths = {};
    for (int v : S) {
        int vPath = 0;
        for (vector<int> u : F[v]) {vPath = vPath + u[k-2];}
        paths[v] = vPath; cntKCol = cntKCol + vPath;
    }
    discrete_distribution<> distr(S.begin(), S.end());
    vector<double> p = distr.probabilities();

    for (int u = 0; u < S.size(); u++) {
        p[u] = paths[S[u]]/cntKCol;
    }
    float successTimes = 0;
    int adj = 0;
    for (int i = 0; i < t; i++) {
        int v = S[distr(gen)];
        vector<node> vertexSubset; 
        vertexSubset.reserve(G[v].outNeighbours.size());

        unordered_map<int, int> labelIndexMap = {};

        for (int j = 0; j < G[v].outNeighbours.size(); j++) {
            node vertexCopy = G[G[v].outNeighbours[j]];
            labelIndexMap[vertexCopy.label] = j;
            vertexSubset.push_back(vertexCopy);
        }

        if (!vertexSubset.empty()) {
            vector<int> R = DPPathSampler(induceSubgraph(vertexSubset), labelIndexMap, k-1); 
            if (R.size() == k-1) {
                R.push_back(v);
                sort(R.begin(), R.end());

                vector<node> C;

                for (int i : R) {
                    C.push_back(G[i]);
                };
                
                successTimes = successTimes + checkClique(C);
            }
            else {adj++;}
        }
        else {
            adj++;
        }
    }
    cout << "cntKCol " << cntKCol << endl;
    return (successTimes/(t-adj)) * cntKCol;
}

int cliqueCounter(vector<vector<int>> partition, vector<node> vertices, int k, int t, int r) {
    int nCliques = 0;
    for (auto v : partition[0]) {
        vector<node> vertexSubset;
        if (vertices[v].outNeighbours.size() >= k-1) {
            for (auto u : vertices[v].outNeighbours) {
                vertexSubset.push_back(vertices[u]);
            }
            nCliques = nCliques + pivoter(induceSubgraph(vertexSubset), k-1, true);
        }
    }
    nCliques = nCliques + estimateClique(vertices, partition[1], k, t, r);

    // code for exact counting of dense partition (gives total exact count)
    // for (auto v : partition[1]) {
    //     vector<node> vertexSubset;
    //     if (vertices[v].outNeighbours.size() >= k-1) {
    //         for (auto u : vertices[v].outNeighbours) {
    //             vertexSubset.push_back(vertices[u]);
    //         }
    //         nCliques = nCliques + pivoter(induceSubgraph(vertexSubset), k-1, true);
    //     }
    // }

    return nCliques;
}

int main() {
    //int k = 3;
    //cout << "Initialised k = " << k << endl;

    vector<node> vertices = vector<node>();
    for (int v = 0; v < n; v++) {
        node curr;
        curr.label = v; curr.deg = 0; curr.colour = n+1;
        vertices.push_back(curr);
    }

    auto t1 = chrono::high_resolution_clock::now();
    getData(vertices);
    cout << "got data" << endl;
    vector<int> ordering = getDegeneracyOrder(vertices);
    cout << "got order" << endl;
    //setDAGNeighbourhoods(vertices, ordering);
    //vector<vector<int>> partition = partitionVertices(k, vertices);
    int r = greedyColouring(vertices, ordering);
    cout << r << "-Coloured \n";
    //vector<int> colourOrdering = getColourOrder(vertices);
    //setColourDAGNeighbourhoods(vertices, colourOrdering);
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Preprocessing Done in (" << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms)" <<endl;

    //auto t15 = chrono::high_resolution_clock::now();
    //vector<vector<int>> H = DPPathCount(vertices, k);
    //auto t16 = chrono::high_resolution_clock::now();
    //cout << "DP Complete (" << chrono::duration_cast<std::chrono::milliseconds>(t16-t15).count() << "ms)" <<endl;
    //DPPathSampling(vertices, k);

    //cout << "Exact " << k << "-clique count: " << pivoter(vertices, k, false) << "\n";
    //cout << cliqueCounter(partition, vertices, k, 10000, r) << "\n";
    testSuffix(vertices, ordering);
}