#include "utils.h"

using namespace std;

void printVector(const vector<int> vec) {
    for (auto i : vec) {
        cout << i << " ";
    }
    cout << "\n";
}

void printMatrix(const vector<vector<int>> matrix) {
    for (auto &row : matrix) {
        for (auto &column : row) {
            cout << column << " ";
        }
        cout << "\n";
    }
}

int choose(int r, int k) {
    if (k == 0) {return 1;}
    return (r * choose(r - 1, k - 1)) / k;
}

vector<int> sortedIntersection(vector<int> v1, vector<int> v2){
    //sort(v1.begin(), v1.end()); sort(v2.begin(), v2.end());
    vector<int> v3;
    set_intersection(v1.begin(), v1.end(), v2.begin(),v2.end(), 
                     back_inserter(v3));
    return v3;
}

vector<int> sortedDifference(vector<int> v1, vector<int> v2){
    //sort(v1.begin(), v1.end()); sort(v2.begin(), v2.end());
    vector<int> v3;
    set_difference(v1.begin(), v1.end(), 
                       v2.begin(), v2.end(), 
                       inserter(v3, v3.end()));
    return v3;
}

vector<int> intersection(vector<int> v1, vector<int> v2){
    sort(v1.begin(), v1.end()); sort(v2.begin(), v2.end());
    vector<int> v3;
    set_intersection(v1.begin(), v1.end(), v2.begin(),v2.end(), 
                     back_inserter(v3));
    return v3;
}

vector<int> difference(vector<int> v1, vector<int> v2){
    sort(v1.begin(), v1.end()); sort(v2.begin(), v2.end());
    vector<int> v3;
    set_difference(v1.begin(), v1.end(), 
                       v2.begin(), v2.end(), 
                       inserter(v3, v3.end()));
    return v3;
}