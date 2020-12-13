#include <vector>
#include <cmath>
#include <vecOps.h>
#include <iostream>
using namespace std;

class solClass {
    public:
        vector<vector<double>> A;
        vector<double> b;
        vector<double> cj;
        vector<const char> x;
        vector<const char> xB;
        vector<double> zj;
        vector<double> zc;
}

solClass simplex(solClass) {
    vector<vector<double>> A = solClass.A;
    vector<double> b = solClass.b;
    vector<double> cj = solClass.cj;
    vector<const char> x = solClass.x;
    vector<const char> xB = solClass.xB;
    vector<double> zj = solClass.zj;
    vector<double> zc = solClass.zc;

    std::cout << xB[0];
    return solution;
}

int main() {
    solClass solution;
    solution.A = {{3, 2, 1, 1, 0, 0,}, {1, 3, 5, 0, 1, 0}, {2, 3, 5, 0, 0, 1}};
    solution.b = {50, 30, 30};
    solution.cj = {10, 5, 15, 0, 0, 0};
    solution.x = {"x1", "x2", "x3", "s1", "s2", "s3"};
    solution.xB = {"s1", "s2", "s3"};
    solution.zj = {0, 0, 0, 0, 0, 0};
    solution.zc = {0, 0, 0, 0, 0, 0};

    solClass solUpdated = simplex(solution);
    return 1;
}