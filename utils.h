#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <tuple>
#include <memory>
#include <algorithm>
#include <queue>
#include <stack>
#include <tuple>
#include <unordered_map>

void printVector(const std::vector<int> vec);

void printMatrix(const std::vector<std::vector<int>> matrix);

int choose(int r, int k);

std::vector<int> sortedIntersection(const std::vector<int> v1, const std::vector<int> v2);

std::vector<int> sortedDifference(const std::vector<int> v1, const std::vector<int> v2);

std::vector<int> intersection(const std::vector<int>, const std::vector<int> v2);

std::vector<int> difference(const std::vector<int>, const std::vector<int> v2);

#endif