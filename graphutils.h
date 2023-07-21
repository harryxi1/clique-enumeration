#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include <fstream>
#include <string>
#include <sstream>

#include "utils.h"

struct node {
int label;
int deg;
int colour; 

std::vector<int> neighbours; // N(v) in undirected graph

std::vector<int> outNeighbours; // N^+(v) in DAG formed by degeneracy order

std::vector<int> colourOutneighbours; // N^+(v) in DAG formed by total colour order
};

std::vector<int> getDegeneracyOrder(std::vector<node> vertices);

std::vector<node> induceSubgraph(std::vector<node> Sv);

float getEdgeDensity(std::vector<node> Sv);

int checkClique(std::vector<node> C);

#endif
