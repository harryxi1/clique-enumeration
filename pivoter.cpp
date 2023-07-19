#include "utils.h"
#include "graphutils.h"
#include "pivoter.h"

using namespace std;

SCT SCTBuilder(vector<node> vertices) {
    SCT T = {vector<int>(), vector<shared_ptr<SCT>>(), vector<tuple<int, bool>>()};
    queue<shared_ptr<SCT>> Q;

    int nVertices = vertices.size();
    T.label.reserve(nVertices); T.children.reserve(nVertices); T.links.reserve(nVertices);

    for (int i = 0; i < nVertices; i++) {
        T.label.push_back(i);

        shared_ptr<SCT> child = make_shared<SCT>(SCT{vertices[i].outneighbours, vector<shared_ptr<SCT>>(), vector<tuple<int, bool>>()});
        T.children.push_back(child);
        T.links.emplace_back(i, 0);
        Q.push(move(child));
    };

    while (!Q.empty()) {
        shared_ptr<SCT> curr = Q.front(); Q.pop();
        if ((*curr).label.empty()) {continue;}
        int p = (*curr).label[0]; int pval = sortedIntersection((*curr).label, vertices[(*curr).label[0]].neighbours).size();

        for (int v = 1; v < (*curr).label.size(); v++) {
            int val = sortedIntersection((*curr).label, vertices[(*curr).label[v]].neighbours).size();
            if (val > pval) {
                pval = val;
                p = (*curr).label[v];
            };
        }
        shared_ptr<SCT> child = make_shared<SCT>(sortedIntersection((*curr).label, vertices[p].neighbours), vector<shared_ptr<SCT>>(), vector<tuple<int, bool>>());
        (*curr).links.emplace_back(p, 1);
        (*curr).children.push_back(child);
        Q.push(move(child));

        vector<int> unionNeighbourhood = vertices[p].neighbours;
        unionNeighbourhood.insert(upper_bound(unionNeighbourhood.begin(), unionNeighbourhood.end(), p), p);
        // sort(unionNeighbourhood.begin(), unionNeighbourhood.end());
        vector<int> diff = sortedDifference((*curr).label, unionNeighbourhood);

        int diffSize = diff.size();
        (*curr).children.reserve(diffSize+1);
        for (int i = 0; i < diffSize; i++) {
            vector<int> neighbourhoodi = sortedIntersection((*curr).label, vertices[diff[i]].neighbours);
            vector<int> sliced(diff.begin(), diff.begin() + i);
            vector<int> diff2 = sortedDifference(neighbourhoodi, sliced);

            shared_ptr<SCT> child = make_shared<SCT> (SCT{diff2, vector<shared_ptr<SCT>>(), vector<tuple<int, bool>>()});
            (*curr).children.push_back(child);
            (*curr).links.emplace_back(diff[i], 0);

            Q.push(move(child));
        }
    }
    return T;
}

void getRootToLeaf(shared_ptr<SCT> root, tuple<int, bool> edgeLink, vector<int> &curr, vector<vector<int>> &total, bool start) {
    if (!start) {
        curr[get<1>(edgeLink)]++;     
    }
    // if current node is leaf
    if ((*root).children.empty()) {
        vector<int> currCopy = curr;
        total.push_back(currCopy);

        curr[get<1>(edgeLink)]--;
        return;
    }
    for (int i = 0; i < (*root).children.size(); i++) {
        vector<int> currCopy = curr;
        getRootToLeaf((*root).children[i], (*root).links[i], curr, total, false);
        curr = currCopy;
    }
}

int pivoter(vector<node> vertices, int k) {
    shared_ptr<SCT> T = make_shared<SCT>(SCTBuilder(vertices));

    int cnt = 0;
    vector<vector<int>> total;
    vector<int> initCurr = {0, 0};
    getRootToLeaf(T, make_tuple(0, 0), initCurr, total, true);
    for (auto T : total) {
        for (int i = 0; i <= T[1]; i++) {
            if (T[0] + i == k) {cnt = cnt + choose(T[1], i);}
        }
    }
    return cnt;
}
