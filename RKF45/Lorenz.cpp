// From matplotlib-cpp-git AUR package
#include <matplotlib-cpp/matplotlibcpp.h>
namespace plt = matplotlibcpp;
#include <iostream>
// Required for setprecision
#include <iomanip>
#include <string>
#include <fstream>
// RKF45 function
#include "RKF45.h"
// Timer
#include <chrono>
using namespace std::chrono;

/**
 * Find dt times the RHS of the ODE expressed as a system of first-order
 * equations.
 *
 * @param params   A vector of parameters.
 * @param t        Time value.
 * @param vars     A vector of dependent variable values.
 * @param dt       Step size.
 * @return         vector of differentials
 */
vector<double> Lorenz(vector<double> params, double t, vector<double> vars,
double dt) {
    // Extract parameters
    double sigma = params[0];
    double beta = params[1];
    double rho = params[2];

    // Extract variable values
    double x = vars[0];
    double y = vars[1];
    double z = vars[2];

    // Derivatives
    double dxdt = sigma*(y-x);
    double dydt = x*(rho-z) - y;
    double dzdt = x*y - beta*z;

    // Differentials
    double dx = dt * dxdt;
    double dy = dt * dydt;
    double dz = dt * dzdt;
    return {dx, dy, dz};
}

/**
 * @brief Solves the problem and provides desired output, such as saved plots
 * and data in a textfile.
 */
int main() {
    // Solution parameters
    double dtInitial = 0.1;
    double epsilon = 1e-9;

    // Initial conditions and domain of integration
    double x0 = 1.0;
    double y0 = 1.0;
    double z0 = 1.0;
    vector<double> conds = {x0, y0, z0};
    double t0 = 0;
    double tf = 200;

    // Problem parameters
    double sigma = 10.0;
    double beta = 8.0/3.0;
    double rho = 28.0;
    vector<double> params = {sigma, beta, rho};

    // Solve problem
    auto start = high_resolution_clock::now(); 
    solClass solution = RKF45(Lorenz, dtInitial, epsilon, params, t0, tf, 
    conds);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "RKF45 took " << duration.count() << " ms to run" << endl;
    vector<double> t = solution.t;
    int k = t.size();
    vector<vector<double>> vars = solution.vars;
    vector<double> x = vars[0];
    vector<double> y = vars[1];
    vector<double> z = vars[2];

    // Write to file
    ofstream myfile;
    myfile.open("Lorenz.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(k), ' ');
    myfile << "t" << std::string(20, ' ');
    myfile << "x" << std::string(17, ' ');
    myfile << "y" << std::string(16, ' ');
    myfile << "z" << "\n";
    // Contents
    for (int i = 0 ; i < k; i++) {
        if (i == 0) {
            myfile << i << std::string(1 + (int)log10(k), ' ');
        } else {
            myfile << i << std::string(1 + (int)log10(k) - (int)log10(i), ' ');
        }
        myfile << setprecision(15) << t[i] << " ";
        myfile << setprecision(15) << x[i] << " ";
        myfile << setprecision(15) << y[i] << "\n";
        myfile << setprecision(15) << z[i] << "\n";
    }

    // Plot using matplotlibcpp
    // You will get linting errors for plt::plot, but no build errors if your
    // matplotlibcpp package is installed and set up properly
    plt::figure(1);
    plt::plot3(x, y, z);
    plt::xlabel("x");
    plt::ylabel("y");
    plt::set_zlabel("z");
    plt::save("Lorenz 3D phaseplot.svg");
    plt::figure(3);
    plt::clf();
    plt::plot(t, x, {{"label", "$x$"}});
    plt::plot(t, y, {{"label", "$y$"}});
    plt::plot(t, z, {{"label", "$z$"}});
    plt::xlabel("Time");
    plt::legend();
    plt::save("Lorenz x, y and z against t plot.svg");
    return 1;
}
