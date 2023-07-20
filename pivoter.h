#ifndef PIVOTER_H
#define PIVOTER_H

#include "utils.h"

struct SCT {
    std::vector<int> label;
    std::vector<std::shared_ptr<SCT>> children;
    std::vector<std::tuple<int, bool>> links;

    SCT(std::vector<int> label, std::vector<std::shared_ptr<SCT>> children, std::vector<std::tuple<int, bool>> links) 
    : label(std::move(label)), children(move(children)), links(move(links)) {}
};

SCT SCTBuilder(std::vector<node> vertices, bool debug);

void getRootToLeaf(std::shared_ptr<SCT> root, std::tuple<int, bool> edgeLink, std::vector<int> &curr, std::vector<std::vector<int>> &total, bool start);

int pivoter(std::vector<node> vertices, int k, bool debug);

#endif