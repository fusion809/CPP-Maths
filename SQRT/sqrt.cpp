// Copyright (C) 2016-2017 Brenton Horne
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <utility>
#include <unistd.h>
#include <algorithm>

using namespace std;

double N    = 100;   // Number of steps in Newton's iteration
int prec    = 15;    // Number of decimal places the sqrt is to be calculated to
double sq   = 3.0;   // Number to square root
double x0   = 1.73;  // Initial x to calculate the root from

int main() {
    std::vector<double> x;
    x.push_back(x0);

    ofstream myfile;
    myfile.open("sqrt.txt");
    myfile.precision(prec);

    myfile << 0;
    myfile << " ";
    myfile << x0;
    for(int i = 1; i<=N; i++)
    {
         x.push_back((1/2.0) * ( x[i-1] + sq/(x[i-1]) ) );
         myfile << "\n";
         myfile << i;
         myfile << " ";
         myfile << x[i];
         usleep(1000);
    }

    myfile.close();
    std::cout.precision(prec);
    std::cout << x[N];
    return 0;
}
