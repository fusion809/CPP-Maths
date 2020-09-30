// This program solves:
// d2y/dx2 + kx dy/dx + x^2 y + x cos(x) = 0 
// Copyright (C) 2017 Brenton Horne
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <utility>
#include <unistd.h>
#include <algorithm>

using namespace std;

double N     = 1000000.0;                 // Number of domain elements we're integrating over
double x0    = 0.0;                       // Starting x
double x1    = 20;                        // End x
double Y0    = 0.0;
double dy0   = 0.0;
double k     = 1;
double h     = (x1 - x0) / double(N);    // step size

double dy2(double k, double x, double y, double dy)
{
    double RHS=(-k*x*dy-x*x*y-x*cos(x));
    return RHS;
}

std::pair<double, double> RK4(double x, double y, double dy, double h)
{
    // kn values are diffs in dy/dx.
    // ln values are diffs in y.
    double k1    = h * dy2(k, x,       y,        dy);
    double l1    = h * dy;

    double k2    = h * dy2(k, x + h/2, y + l1/2, dy + k1/2);
    double l2    = h * (dy + k1/2);

    double k3    = h * dy2(k, x + h/2, y + l2/2, dy + k2/2);
    double l3    = h * (dy + k2/2);

    double k4    = h * dy2(k, x + h,   y + l3,   dy + k3);
    double l4    = h * (dy + k3);

    double diff1 = (l1 + 2 * l2 + 2 * l3 + l4)/float(6);  // diff in x.
    double diff2 = (k1 + 2 * k2 + 2 * k3 + k4)/float(6);  // diff in y.

    return {diff1, diff2};
}

int main()
{
    std::vector<double> x;
    x.push_back(x0);

    std::vector<double> y;
    y.push_back(Y0);

    std::vector<double> dy;
    dy.push_back(dy0);

    ofstream myfile;
    myfile.open("Sol-y.txt");

    ofstream myfiledy;
    myfiledy.open("Sol-dy.txt");

    ofstream myfilephase;
    myfilephase.open("Sol-phase.txt");

    double miny = y[0];

    for(int i = 1; i<=N; i++)
    {
        std::pair<double, double> diff = RK4(x[i-1],y[i-1],dy[i-1],h);

        x.push_back(   x[i-1] + h);
        y.push_back(   y[i-1] + diff.first  );
        dy.push_back( dy[i-1] + diff.second );
        
        myfile      << x[i-1];
        myfile      << " "     << y[i-1]   << "\n";

        myfiledy    << x[i-1];
        myfiledy    << " "     << dy[i-1]  << "\n";

        myfilephase << y[i-1];
        myfilephase << " "     << dy[i-1]  << "\n";

        usleep(1000);
    }

    myfile   << x[N];
    myfile   << " "   << y[N];

    myfiledy << x[N];
    myfiledy << " "   << dy[N];

    myfile.close();
    myfiledy.close();

    system("gnuplot -p RandomExample2.gp");
    system("rsvg-convert -w 2000 -o RandomExample2.png RandomExample2.svg");
    system("rsvg-convert -w 2000 -o RandomExample2-dy.png RandomExample2-dy.svg");
    system("rsvg-convert -w 2000 -o RandomExample2-phase.png RandomExample2-phase.svg");

    return 0;
}
