#include "graphutils.h"

std::vector<int> getDegeneracyOrder(std::vector<node> vertices) {
    std::vector<int> L = std::vector<int>();
    std::vector<int> dv = std::vector<int>();
    std::vector<std::vector<int>> D;
    int k = 0;

    int max_deg = 0;
    for (int v = 0; v < vertices.size(); v++) {
        int curr_deg = vertices[v].deg;
        dv.push_back(curr_deg);
        if (curr_deg > max_deg) {max_deg = curr_deg;}
    }
    for (int i = 0; i < max_deg+1; i++) {
        D.push_back({});
    }
    for (int v = 0; v < vertices.size(); v++) {
        D[vertices[v].deg].push_back(v);
    }
    for (int t = 0; t < vertices.size(); t++) {
        int i = 0;
        while (D[i].empty()) {i++;}
        k = std::max(k, i);
        //L.insert(L.begin(), D[i][0]);
        L.push_back(D[i][0]);
        D[i].erase(D[i].begin());
        for (int j = 0; j < vertices[L[0]].neighbours.size(); j++) {
            int curr_neigh = vertices[L[0]].neighbours[j];
            int neigh_dv = dv[curr_neigh];
            if (find(L.begin(), L.end(), curr_neigh) == L.end()) {
                D[neigh_dv].erase(remove(D[neigh_dv].begin(), D[neigh_dv].end(), curr_neigh), D[neigh_dv].end());
                dv[curr_neigh]--;
                D[dv[curr_neigh]].push_back(curr_neigh);
                sort(D[dv[curr_neigh]].begin(), D[dv[curr_neigh]].end());
            }
        }
    }
    return L;
}

std::vector<node> induceSubgraph(std::vector<node> Sv) {
    // given set of vertices
    int size = Sv.size();
    std::vector<int> SvLabel; SvLabel.reserve(size);
    for (node v : Sv) {
        SvLabel.push_back(v.label);
    }
    std::vector<node> inducedSubgraph; inducedSubgraph.reserve(size);
    for (node v : Sv) {
        node vertexCopy = v;
        vertexCopy.neighbours = sortedIntersection(SvLabel, vertexCopy.neighbours);
        vertexCopy.outNeighbours = sortedIntersection(SvLabel, vertexCopy.outNeighbours);
        vertexCopy.colourOutneighbours = sortedIntersection(SvLabel, vertexCopy.colourOutneighbours);
        inducedSubgraph.push_back(vertexCopy);
    }
    return inducedSubgraph;
}

float getEdgeDensity(std::vector<node> Sv) {
    std::vector<node> inducedSubgraph = induceSubgraph(Sv);
    float degsum = 0;
    for (node v : inducedSubgraph) {degsum = degsum + v.neighbours.size();}
    return (degsum/2)/(choose(Sv.size(), 2));
}

int checkClique(std::vector<node> C) {
    return getEdgeDensity(C) == 1;
}