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

std::vector<int> inneighbours; // N^-(v) in DAG formed by degeneracy order
std::vector<int> outneighbours; // N^+(v) in DAG formed by degeneracy order

std::vector<int> colour_inneighbours; // N^-(v) in DAG formed by total colour order
std::vector<int> colour_outneighbours; // N^+(v) in DAG formed by total colour order
};

#endif
