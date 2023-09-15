#include "stdafx.h"

#include <vector>
#include <numeric>
#include <chrono>

#include "plf_reorderase.h"
#include "hash_table6.hpp"
#include "Snap.h"

// data directories
#define FILENAME "testData/com-youtube.ungraph.txt"

// precomputed graph properties
#define DEGENERACY 51

using namespace std;

float choose(const int r, const int k) {
    if (k == 0) { return 1; }
    return (r * choose(r - 1, k - 1)) / k;
}

emhash6::HashMap<int, float> getCombValues(const int k) {
    emhash6::HashMap<int, float> combValues;
    for (int i = 0; i <= DEGENERACY + 1; ++i) {
        combValues[i] = choose(i, k);
    }
    return combValues;
}

int getMaxDegree(const PUNGraph G) {
    return (G->GetNI(TSnap::GetMxDegNId(G))).GetDeg();
}

struct degeneracyObject {
    vector<int> degOrder;
    emhash6::HashMap<int, int> degIndexMap;
    int degeneracy;

    // constructor
    explicit degeneracyObject(const int n) : degOrder{ vector<int>(n, 0) },
                                  degIndexMap{ emhash6::HashMap<int, int>() }, 
                                  degeneracy{ 0 } {}
};

degeneracyObject computeDegeneracyOrder(const PUNGraph G) {
    const int n = G->GetNodes();
    degeneracyObject output{n};

    vector<vector<int>> D(getMaxDegree(G) + 1, vector<int>()); // buckets
    emhash6::HashMap<int, int> dv; // node id to bucket
    emhash6::HashMap<int, int> nodeLocator; // node id to index (in D[i])

    int minBucket = G->BegNI().GetDeg();

    // initialise containers - traverse all nodes in G
    int currID; int currDeg;
    for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        currID = NI.GetId(); currDeg = NI.GetDeg();

        dv[currID] = currDeg;
        D[currDeg].push_back(currID);

        // nodeLocator[currID] = D[currDeg].size() - 1;
        nodeLocator.insert_unique(currID, D[currDeg].size() - 1);
        minBucket = min(minBucket, currDeg);
    }

    int neighbourID; // variable for holding current neighbourid

    // repeat n times
    for (int i = 0; i < n; ++i) {
        while (D[minBucket].empty()) { ++minBucket; }
        output.degeneracy = max(output.degeneracy, minBucket);

        // select last vertex in bucket, place in order, then remove
        output.degOrder[i] = (D[minBucket].back());
        D[minBucket].pop_back(); dv[output.degOrder[i]] = -1;
        output.degIndexMap[output.degOrder[i]] = i;
        minBucket = max(minBucket - 1, 0);

        // reduce neighbours buckets
        TUNGraph::TNodeI NI = G->GetNI(output.degOrder[i]);
        for (int j = 0; j < NI.GetDeg(); ++j) {
            if (dv[NI.GetNbrNId(j)] != -1) {
                neighbourID = NI.GetNbrNId(j);

                // remove from buckets
                if (D[dv[neighbourID]].back() == neighbourID) {
                    D[dv[neighbourID]].pop_back();
                }
                else {
                    nodeLocator[D[dv[neighbourID]].back()] = nodeLocator[neighbourID];

                    plf::reorderase(D[dv[neighbourID]],
                        D[dv[neighbourID]].begin() + nodeLocator[neighbourID]);
                }
                --dv[neighbourID];
                D[dv[neighbourID]].push_back(neighbourID);
                nodeLocator[neighbourID] = D[dv[neighbourID]].size() - 1;
            }
        }
    }
    return output;
}

TIntV getOutNeighbourhoodID(
    emhash6::HashMap<int, int>& degIndexMap, const TUNGraph::TNodeI NI) {
    TIntV outNeighbourhood; 
    for (int i = 0; i < NI.GetDeg(); ++i) {
        if (degIndexMap[(NI.GetId())] < degIndexMap[(NI.GetNbrNId(i))]) {
            outNeighbourhood + NI.GetNbrNId(i);
        }
    }
    return outNeighbourhood;
}

bool checkInClique(int nodeID, vector<TUNGraph::TNodeI>& currClique) {
    for (int i = currClique.size() - 1; i >= 0; --i) {
        if (!currClique[i].IsNbrNId(nodeID)) {return 0;}
    }
    return 1;
}

float countCliqueSuffix(const PUNGraph G, const int k, const float mu) {
    if (k == 2) {return G->GetEdges();}
    // TODO: calculate this once
    emhash6::HashMap<int, float> combValues = getCombValues(k-1);

    degeneracyObject degObj = computeDegeneracyOrder(G);
    float total = 0;

    vector<TUNGraph::TNodeI> clique;
    for (TUNGraph::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
        const PUNGraph sG = TSnap::GetSubGraph(G, getOutNeighbourhoodID(degObj.degIndexMap, NI));
        if (!sG->Empty()) {
            // handshaking lemma to get average degree
            int misses = 0;
            if ((2 * (sG->GetEdges()) / (sG->GetNodes())) <= k) {
                total += countCliqueSuffix(sG, k-1, mu);
            }
            else {
                // dense
                degeneracyObject subDegObj = computeDegeneracyOrder(sG);
                for (int j = sG->GetNodes() - 1; j >= 0; --j) {
                    ++misses;
                    if (clique.empty() || checkInClique((subDegObj.degOrder)[j], clique)) {
                        clique.push_back(
                            sG->GetNI((subDegObj.degOrder)[j])
                        );
                        misses = 0;
                    }
                    if (misses > max(mu * sG->GetNodes(), 1)) {
                        break;
                    }
                }
                total += combValues[clique.size()];
                clique.clear();
            }
        }
    }
    return total;
}

int main() {
    const int k = 10;
    printf("Loading Graph:\n");
    PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(FILENAME, 0, 1);

    // if we want to remove zero degree nodes
    TSnap::DelZeroDegNodes(G);

    printf("Counting:\n");

    auto start = chrono::high_resolution_clock::now();
    float total = countCliqueSuffix(G, k, 1);
    auto stop = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::microseconds>(stop - start);
    printf("Completed in %lld microseconds\n", duration.count());
    printf("Found %f, %d-cliques", total, k);

    return 0;
}
