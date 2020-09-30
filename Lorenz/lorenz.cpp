// Lorenz equation solver
// Copyright (C) 2017 fusion809
// Requires C++11 support, as tuple is a type found only in C++11 or later
// Source (of tuples only existing in C++11 and later): http://en.cppreference.com/w/cpp/utility/tuple
// Clang 3.3 and later, as well as GCC 4.7 and later have at least experimental support for C++11
// I have compiled this program successfully with GCC 4.8, 5.4, 6.1, 6.2, 6.3 and 7.2.
#include <tuple>                                           // Needed for tuples, supported in C++11 and later
#include <vector>                                          // Needed for vectors
#include <fstream>                                         // Needed for ofstream

// Standard template library
using namespace std;

// Number of domain elements we're integrating over.
// Cannot go beyond 10,000,000 as otherwise the SVG becomes too large for rsvg-convert to handle.
double N     = 10000.0;

// Integration interval
double t0    = 0.0;                                        // Starting time
double t1    = 100;                                        // End time

// Initial conditions
double X0    = 1.0;                                        // x(t0)
double Y0    = 1.0;                                        // y(t0)
double Z0    = 1.0;                                        // z(t0)

// Parameters, a copy of the Lorenz equations (using these exact symbols to denote the parameters) is in Lorenz_Equations.png
double rho   = 28;
double sigma = 10;
double beta  = 8/3.0;

// Step size
double h     = (t1 - t0) / double(N);

// dX is a function that returns a 3-tuple; first entry is dx/dt, second is dy/dt and third is dz/dt.
std::tuple<double, double, double> dX(double h, double sigma, double rho, double beta, double t, double x, double y, double z)
{
    // Used to use:
    // return {h * sigma * (y - x), h * x * (rho - z) - h * y, h * x * y - h * beta * z}
    // But that is inconsistent with GCC < 6 and Clang.
    // Eugene Porotov (https://www.facebook.com/eugene.porotov?fref=ufi) came up with the idea to use make_tuple instead
    return std::make_tuple(h * sigma * (y - x), h * x * (rho - z) - h * y, h * x * y - h * beta * z);
}

// RK4 is a function that returns a 3-tuple of Runge-Kutta fourth-order differences.
// First entry is the diff in x, second is diff in y and third is diff in z.
std::tuple<double, double, double> RK4(double t, double x, double y, double z, double h)
{
    // kn values are diffs in x.
    // ln values are diffs in y.
    // mn values are diffs in z.
    // diff1
    std::tuple<double, double, double> diff1 = dX(h, sigma, rho, beta, t,       x,        y,       z);
    double k1    = std::get<0>(diff1);
    double l1    = std::get<1>(diff1);
    double m1    = std::get<2>(diff1);
    
    // diff2
    std::tuple<double, double, double> diff2 = dX(h, sigma, rho, beta, t + h/2, x + k1/2, y + l1/2, z + m1/2);
    double k2    = std::get<0>(diff2);
    double l2    = std::get<1>(diff2);
    double m2    = std::get<2>(diff2);

    // diff3
    std::tuple<double, double, double> diff3 = dX(h, sigma, rho, beta, t + h/2, x + k2/2, y + l2/2, z + m2/2);
    double k3    = std::get<0>(diff3);
    double l3    = std::get<1>(diff3);
    double m3    = std::get<2>(diff3);

    // diff4
    std::tuple<double, double, double> diff4 = dX(h, sigma, rho, beta, t + h,   x + k3,   y + l3,   z + m3);
    double k4    = std::get<0>(diff4);
    double l4    = std::get<1>(diff4);
    double m4    = std::get<2>(diff4);

    // diff in x, y, z
    double diffx = (k1 + 2 * k2 + 2 * k3 + k4)/float(6);   // diff in x.
    double diffy = (l1 + 2 * l2 + 2 * l3 + l4)/float(6);   // diff in y.
    double diffz = (m1 + 2 * m2 + 2 * m3 + m4)/float(6);   // diff in z.

    return std::make_tuple(diffx, diffy, diffz);
}

int main()
{
    // t vector
    std::vector<double> t;
    // add t0 to vector
    t.push_back(t0);

    // initialize x vector
    std::vector<double> x;
    // add X0 to vector
    x.push_back(X0);

    // initialize y vector
    std::vector<double> y;
    // add Y0 to vector
    y.push_back(Y0);

    // initialize z vector
    std::vector<double> z;
    // add Z0 to vector
    z.push_back(Z0);

    // open filetx ofstream
    ofstream filetx;
    filetx.open("lorenz-tx.txt");

    // open filety ofstream
    ofstream filety;
    filety.open("lorenz-ty.txt");

    // open filetz ofstream
    ofstream filetz;
    filetz.open("lorenz-tz.txt");

    // open filexy ofstream
    ofstream filexy;
    filexy.open("lorenz-xy.txt");

    // open filexz ofstream
    ofstream filexz;
    filexz.open("lorenz-xz.txt");

    // open fileyz ofstream
    ofstream fileyz;
    fileyz.open("lorenz-yz.txt");

    // open filexyz ofstream
    ofstream filexyz;
    filexyz.open("lorenz-xyz.txt");

    // Loop over elements in t
    for(int i = 1; i<=N; i++)
    {
        // delta X
        // Where X = (x, y, z)
        // So delta X is the change in X elements
        std::tuple<double, double, double> deltaX = RK4(t[i-1], x[i-1], y[i-1], z[i-1], h);

        // diffs in x, y and z
        double deltax = std::get<0>(deltaX);
        double deltay = std::get<1>(deltaX);
        double deltaz = std::get<2>(deltaX);

        // add an ith element to t
        t.push_back( t[i-1] + h      );
        // add an ith element to x
        x.push_back( x[i-1] + deltax );
        // add an ith element to y
        y.push_back( y[i-1] + deltay );
        // add an ith element to z
        z.push_back( z[i-1] + deltaz );

        // filetx (i-1)th entry
        filetx    << t[i-1];
        filetx    << " "     << x[i-1]   << "\n";
        filetx.precision(15);

        // filety (i-1)th entry
        filety    << t[i-1];
        filety    << " "     << y[i-1]  << "\n";
        filety.precision(15);

        // filetz (i-1)th entry
        filetz    << t[i-1];
        filetz    << " "     << z[i-1]  << "\n";
        filetz.precision(15);

        // filexy (i-1)th entry
        filexy    << x[i-1];
        filexy    << " "     << y[i-1]  << "\n";
        filexy.precision(15);

        // filexz (i-1)th entry
        filexz    << x[i-1];
        filexz    << " "     << z[i-1]  << "\n";
        filexz.precision(15);

        // fileyz (i-1)th entry
        fileyz    << y[i-1];
        fileyz    << " "     << z[i-1]  << "\n";
        fileyz.precision(15);

        // filexyz (i-1)th entry
        filexyz << x[i-1];
        filexyz << " "     << y[i-1];
        filexyz << " "     << z[i-1]  << "\n";
        filexyz.precision(15);
    }

    // filetx Nth entry
    filetx  << t[N];
    filetx  << " "   << x[N];

    // filety Nth entry
    filety  << t[N];
    filety  << " "   << y[N];

    // filetz Nth entry
    filetz  << t[N];
    filetz  << " "   << z[N];

    // filexy Nth entry
    filexy << x[N];
    filexy << " "   << y[N];

    // filexz Nth entry
    filexz << x[N];
    filexz << " "   << z[N];

    // fileyz Nth entry
    fileyz << y[N];
    fileyz << " "   << z[N];

    // filexyz Nth entry
    filexyz << x[N];
    filexyz << " " << y[N];
    filexyz << " " << z[N];

    // close myfile streams
    filetx.close();
    filety.close();
    filetz.close();
    filexy.close();
    filexz.close();
    fileyz.close();
    filexyz.close();

    // Generate svg plots
    system("gnuplot -p lorenz.gp");

    // convert svg plots to png format
    system("rsvg-convert -w 2000 -o lorenz-tx.png lorenz-tx.svg");
    system("rsvg-convert -w 2000 -o lorenz-ty.png lorenz-ty.svg");
    system("rsvg-convert -w 2000 -o lorenz-tz.png lorenz-tz.svg");
    system("rsvg-convert -w 2000 -o lorenz-xy.png lorenz-xy.svg");
    system("rsvg-convert -w 2000 -o lorenz-xz.png lorenz-xz.svg");
    system("rsvg-convert -w 2000 -o lorenz-yz.png lorenz-yz.svg");
    system("rsvg-convert -w 2000 -o lorenz-xyz.png lorenz-xyz.svg");

    // exit
    return 0;
}
