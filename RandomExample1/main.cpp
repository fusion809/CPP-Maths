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

double N     = 100000.0;                 // Number of domain elements we're integrating over
double t0    = 0.01;                     // Starting time
double t1    = 10;                       // End time
double x0    = 0.0;                      // Initial x value
double dx0   = 1.0;                      // Initial dx value
double k     = 2.0;                      // parameters
double l     = 0.5;
double m     = 1.0;
double n     = 1.0;
double h     = (t1 - t0) / double(N);    // step size

double dx2(double k, double l, double m, double n, double t, double x, double dx)
{
    return - (k/t)*dx - (l/(t*t))*x + m*(t+n)/(t*t);
}

std::pair<double, double> RK4(double t, double x, double dx, double h)
{
    // kn values are diffs in dy/dx.
    // ln values are diffs in x.
    double k1    = h * dx2(k, l, m, n, t,       x,        dx);
    double l1    = h * dx;
    double k2    = h * dx2(k, l, m, n, t + h/2, x + l1/2, dx + k1/2);
    double l2    = h * (dx + k1/2);
    double k3    = h * dx2(k, l, m, n, t + h/2, x + l2/2, dx + k2/2);
    double l3    = h * (dx + k2/2);
    double k4    = h * dx2(k, l, m, n, t + h,   x + l3,   dx + k3);
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
    myfile.open("randomexample1.txt");

    ofstream myfiledx;
    myfiledx.open("randomexample1-dx.txt");

    ofstream myfilephase;
    myfilephase.open("randomexample1-phase.txt");

    double minx = x[0];
    double maxx = x[0];

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
        {
           minx=x[i];
        } else if(x[i]>maxx) 
        {
           maxx=x[i];
        }

        usleep(1000);
    }

    myfile   << t[N];
    myfile   << " "   << x[N];

    myfiledx << t[N];
    myfiledx << " "   << dx[N];

    myfile.close();
    myfiledx.close();

    system("gnuplot -p main.gp");
    system("rsvg-convert -o randomexample1-standard.png -w 2000 randomexample1-standard.svg");
    system("rsvg-convert -o randomexample1-derivatve.png -w 2000 randomexample1-derivatve.svg");
    system("rsvg-convert -o randomexample1-phase.png -w 2000 randomexample1-phase.svg");

    std::cout << "Minimum (x):    " << std::scientific;
    std::cout.precision(15);
    std::cout << minx << "\n";
    std::cout << "Maximum (x):    ";
    std::cout << maxx;

    return 0;
}
