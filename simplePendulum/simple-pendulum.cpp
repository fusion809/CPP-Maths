// Simple Pendulum solver
// This program solves the problem of the simple pendulum, operating under
// Newtonian gravity with the acceleration due to gravity being relatively
// constant at 9.8 m/s^2, with x being the angle from the x axis.
// Copyright (C) 2016-2018 Brenton Horne
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <utility>
#include <unistd.h>
#include <algorithm>

using namespace std;

double N     = 100000.0;                 // Number of domain elements we're integrating over
double t0    = 0.0;                      // Starting time
double t1    = 20;                       // End time
double x0    = 0.0;                      // Initial angle from the x axis
double dx0   = 0.0;                      // Initial rate of change in angle from x axis
double g     = 9.8;                      // acceleration due to gravity in m/s^2
double l     = 1.0;                      // Length of pendulum bob
double h     = (t1 - t0) / double(N);    // step size

double dx2(double g, double l, double t, double x, double dx)
{
    return - ( g / l ) * cos(x);
}

std::pair<double, double> RK4(double t, double x, double dx, double h)
{
    // kn values are diffs in dy/dx.
    // ln values are diffs in x.
    double k1    = h * dx2(g, l, t,       x,        dx);
    double l1    = h * dx;

    double k2    = h * dx2(g, l, t + h/2, x + l1/2, dx + k1/2);
    double l2    = h * (dx + k1/2);

    double k3    = h * dx2(g, l, t + h/2, x + l2/2, dx + k2/2);
    double l3    = h * (dx + k2/2);

    double k4    = h * dx2(g, l, t + h,   x + l3,   dx + k3);
    double l4    = h * (dx + k3);

    double diff1 = (l1 + 2 * l2 + 2 * l3 + l4)/float(6);  // diff in x.
    double diff2 = (k1 + 2 * k2 + 2 * k3 + k4)/float(6);  // diff in y.

    return {diff1, diff2};
}

int main()
{
    std::vector<double> t;
    t.push_back(t0);

    std::vector<double> x;
    x.push_back(x0);

    std::vector<double> dx;
    dx.push_back(dx0);

    ofstream myfile;
    myfile.open("simple-pendulum.txt");

    ofstream myfiledx;
    myfiledx.open("simple-pendulum-dx.txt");

    ofstream myfilephase;
    myfilephase.open("simple-pendulum-phase.txt");

    double minx = x[0];

    for(int i = 1; i<=N; i++)
    {
        std::pair<double, double> diff = RK4(t[i-1],x[i-1],dx[i-1],h);

        t.push_back(   t[i-1] + h);
        x.push_back(   x[i-1] + diff.first  );
        dx.push_back( dx[i-1] + diff.second );

        myfile      << t[i-1];
        myfile      << " "     << x[i-1]   << "\n";

        myfiledx    << t[i-1];
        myfiledx    << " "     << dx[i-1]  << "\n";

        myfilephase << x[i-1];
        myfilephase << " "     << dx[i-1]  << "\n";

        if(x[i]<minx)
        minx=x[i];

        usleep(1000);
    }

    myfile   << t[N];
    myfile   << " "   << x[N];

    myfiledx << t[N];
    myfiledx << " "   << dx[N];

    myfile.close();
    myfiledx.close();

    system("gnuplot -p simple-pendulum.gp");
    system("rsvg-convert -w 2000 -o simple-pendulum-standard.png simple-pendulum-standard.svg");
    system("rsvg-convert -w 2000 -o simple-pendulum-derivative.png simple-pendulum-derivative.svg");
    system("rsvg-convert -w 2000 -o simple-pendulum-phase.png simple-pendulum-phase.svg");

    std::cout << "Error:    " << std::scientific;
    std::cout.precision(15);
    std::cout << minx + M_PI;
    std::cout << ".\nx[N] is: " << x[N] << ".\n";

    return 0;
}
