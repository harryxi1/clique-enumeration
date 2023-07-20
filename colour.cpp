#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
#include <random>
#include <chrono>

#include "utils.h"
#include "graphutils.h"
#include "pivoter.h"

using namespace std;

random_device rd;
mt19937 gen(rd());

#define FILENAME "data/email_EU_edges.csv"
#define n 1005
#define m 25571

vector<vector<int>> getData(vector<node> &vertices){
    int count = 0;

    string filename = FILENAME;
    string line;

    ifstream data;
    data.open(filename);
    vector<vector<int>> adjacencyMatrix(n, vector<int>(n, 0));

    while (getline(data, line) && (count < m + 1)) {
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

            adjacencyMatrix[u][v] = 1;
            adjacencyMatrix[v][u] = 1;
        }
    }
    return adjacencyMatrix;
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
        vertices[v].neighbours.erase(unique(vertices[v].neighbours.begin(), vertices[v].neighbours.end() ), vertices[v].neighbours.end());
        vertices[v].outNeighbours.erase(unique(vertices[v].outNeighbours.begin(), vertices[v].outNeighbours.end() ), vertices[v].outNeighbours.end());
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

vector<vector<int>> DPPathCount(vector<node> vertices, int k) {
    vector<vector<int>> H(vertices.size(), vector<int>(k, 0));
    for (int i = 0; i < vertices.size(); i++) {H[i][0] = 1;}
    for (int j = 1; j < k; j++) {
        for (int i = 0; i < vertices.size(); i++) {
            for (int neigh = 0; neigh < vertices[i].colourOutneighbours.size(); neigh++) {
                node currNeighbour = vertices[vertices[i].colourOutneighbours[neigh]];
                H[i][j] = H[i][j] + H[currNeighbour.label][j-1];
            }
        }
    }
    return H;
}

vector<int> DPPathSampling(vector<node> vertices, int k) {
    vector<vector<int>> H = DPPathCount(vertices, k);
    vector<int> R = {};
    vector<int> Q = {};
    for (int i = 0; i < vertices.size(); i++) {Q.push_back(vertices[i].label);}
    for (int i = 0; i < k-1; i++) {
        int cnt = 0;
        for (int u = 0; u < Q.size(); u++) {cnt = cnt + H[Q[u]][k - i - 1];}
        cout << cnt << endl;
        discrete_distribution<> distr(Q.begin(), Q.end());
        vector<double> p = distr.probabilities();
        for (int u = 0; u < Q.size(); u++) {p[u] = H[Q[u]][k - i - 1]/cnt;}
        int u = distr(gen);
        R.push_back(u); Q = vertices[u].colourOutneighbours;
    }
    return R;
}

float estimateClique(vector<node> vertices, vector<int> S, int k, int t, int r) {
    unordered_map<int, vector<vector<int>>> F = {};
    for (auto v : S) {
        // build induced subgraph w/o function
        vector<node> inducedSubgraph;
        inducedSubgraph.reserve(vertices[v].outNeighbours.size());
        for (auto u : vertices[v].outNeighbours) {
            node vertexCopy = vertices[u];
            vertexCopy.colourOutneighbours = intersection(vertices[v].outNeighbours , vertexCopy.colourOutneighbours);
            inducedSubgraph.push_back(vertexCopy);
        }
        F[v] = DPPathCount(inducedSubgraph, k-1);
    }
    int cntKCol = 0;
    for (auto v : S) {cntKCol = cntKCol + F[v][r-1][k-2];}
    
    discrete_distribution<> distr(S.begin(), S.end());
    vector<double> p = distr.probabilities();
    for (int u = 0; u < S.size(); u++) {p[S[u]] = F[S[u]][r-1][k-2]/cntKCol;}
    int successTimes = 0;
    for (int i = 0; i < t; i++) {
        int v = distr(gen);
        // build induced subgraph w/o function
        vector<node> inducedSubgraph;
        inducedSubgraph.reserve(vertices[v].outNeighbours.size());
        for (auto u : vertices[v].outNeighbours) {
            node vertexCopy = vertices[u];
            vertexCopy.colourOutneighbours = intersection(vertices[v].outNeighbours , vertexCopy.colourOutneighbours);
            inducedSubgraph.push_back(vertexCopy);
        }

        vector<int> R = DPPathSampling(inducedSubgraph, k-1); R.push_back(v);
        vector<node> C;
        for (int i : R) {
            C.push_back(vertices[i]);
        }; successTimes = successTimes + checkClique(C);
    }
    float pk = successTimes/t;
    return pk * cntKCol;
}

int checkSuffix(vector<node> vertices, vector<int> ordering) {
    vector<node> suffix;
    for (int i = ordering.size()-1; i <= 0; i--) {
        suffix.push_back(vertices[i]);
        if (!checkClique(suffix)) {break;}
    }
    return suffix.size();
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

    // for (auto v : partition[1]) {
    //     vector<node> vertexSubset;
    //     if (vertices[v].outNeighbours.size() >= k-1) {
    //         for (auto u : vertices[v].outNeighbours) {
    //             vertexSubset.push_back(vertices[u]);
    //         }
    //         nCliques = nCliques + pivoter(induceSubgraph(vertexSubset), k-1, true);
    //     }
    // }
    nCliques = nCliques + estimateClique(vertices, partition[1], k, 100000, r);

    return nCliques;
}

int main() {
    int k = 3;
    cout << "Initialised k = " << k << endl;

    vector<node> vertices = vector<node>();
    for (int v = 0; v < n; v++) {
        node curr;
        curr.label = v; curr.deg = 0; curr.colour = n+1;
        vertices.push_back(curr);
    }

    auto t1 = chrono::high_resolution_clock::now();
    vector<vector<int>> adjacencyMatrix = getData(vertices);
    auto t2 = chrono::high_resolution_clock::now();
    cout << "Data Retrieved (" << chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << "ms)" <<endl;

    auto t3 = chrono::high_resolution_clock::now();
    vector<int> ordering = getDegeneracyOrder(vertices);
    auto t4 = chrono::high_resolution_clock::now();
    cout << "Vertices Ordered (" << chrono::duration_cast<std::chrono::milliseconds>(t4-t3).count() << "ms)" <<endl;

    auto t5 = chrono::high_resolution_clock::now();
    setDAGNeighbourhoods(vertices, ordering);
    auto t6 = chrono::high_resolution_clock::now();
    cout << "Formed DAG (" << chrono::duration_cast<std::chrono::milliseconds>(t6-t5).count() << "ms)" <<endl;

    auto t7 = chrono::high_resolution_clock::now();
    vector<vector<int>> partition = partitionVertices(k, vertices);
    auto t8 = chrono::high_resolution_clock::now();
    cout << "V Partitioned (" << chrono::duration_cast<std::chrono::milliseconds>(t8-t7).count() << "ms)" <<endl;
    cout << partition[0].size() << "\n";
    auto t9 = chrono::high_resolution_clock::now();
    int r = greedyColouring(vertices, ordering);
    auto t10 = chrono::high_resolution_clock::now();
    cout << r << "-Coloured (" << chrono::duration_cast<std::chrono::milliseconds>(t10-t9).count() << "ms)" <<endl;

    auto t11 = chrono::high_resolution_clock::now();
    vector<int> colourOrdering = getColourOrder(vertices);
    auto t12 = chrono::high_resolution_clock::now();
    cout << "Colour Ordered (" << chrono::duration_cast<std::chrono::milliseconds>(t12-t11).count() << "ms)" <<endl;

    auto t13 = chrono::high_resolution_clock::now();
    setColourDAGNeighbourhoods(vertices, colourOrdering);
    auto t14 = chrono::high_resolution_clock::now();
    cout << "Colour DAG Formed (" << chrono::duration_cast<std::chrono::milliseconds>(t14-t13).count() << "ms)" <<endl;

    auto t15 = chrono::high_resolution_clock::now();
    vector<vector<int>> H = DPPathCount(vertices, k);
    auto t16 = chrono::high_resolution_clock::now();
    cout << "DP Complete (" << chrono::duration_cast<std::chrono::milliseconds>(t16-t15).count() << "ms)" <<endl;
    DPPathSampling(vertices, k);

    cout << pivoter(vertices, k, false) << "\n";
    cout << cliqueCounter(partition, vertices, k, 10, r) << "\n";
}