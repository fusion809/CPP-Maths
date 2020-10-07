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

/**
 * Find dt times the RHS of the ODE expressed as a system of first-order 
 * equations.
 *
 * @param params   A vector of parameters.
 * @param t        Time value.
 * @param vars     A vector of dependent variable values.
 * @param dt       Step size.
 * @return         vector of dtheta, dthetaDot
 */
vector<double> simpPen(vector<double> params, double t, vector<double> vars, 
double dt) {
    double g = params[0];
    double l = params[1];
    double theta = vars[0];
    double thetaDot = vars[1];
    double thetaDDot = -g/l * cos(theta);
    return {dt*thetaDot, dt*thetaDDot};
}

/**
 * @brief Solves the problem and provides desired output, such as saved plots 
 * and data in a textfile.
 */
int main() {
    // Solution parameters
    double epsilon = 1e-11;
    double dtInitial = 0.1;

    // Initial conditions and domain of integration
    double theta0 = 0;
    double thetaDot0 = 0;
    vector<double> conds = {theta0, thetaDot0};
    double t0 = 0;
    double tf = 10;

    // Problem parameters
    double g = 9.81;
    double l = 1.0;
    vector<double> params = {g, l};

    // Solve problem
    solClass solution = RKF45(simpPen, dtInitial, epsilon, params, t0, tf, 
    conds);
    vector<double> t = solution.t;
    int k = t.size();
    vector<vector<double>> vars = solution.vars;
    vector<double> theta = vars[0];
    vector<double> thetaDot = vars[1];

    // Write to file
    ofstream myfile;
    myfile.open("simplePendulum.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(k), ' ');
    myfile << "t" << std::string(19, ' ');
    myfile << "theta" << std::string(17, ' ');
    myfile << "thetaDot" << "\n";
    // Contents
    for (int i = 0 ; i < k; i++) {
        if (i == 0) {
            myfile << i << std::string(1 + (int)log10(k), ' ');
        } else {
            myfile << i << std::string(1 + (int)log10(k) - (int)log10(i), ' ');
        }
        myfile << setprecision(15) << t[i] << " ";
        myfile << setprecision(15) << theta[i] << " ";
        myfile << setprecision(15) << thetaDot[i] << "\n";
    }

    // Plot using matplotlibcpp
    // You will get linting errors for plt::plot, but no build errors if your
    // matplotlibcpp package is installed and set up properly
    plt::figure(1);
    plt::plot(t, theta, {{"label", "$\\theta$"}});
    plt::plot(t, thetaDot, {{"label", "$\\dot{\\theta}$"}});
    plt::xlabel("$t$");
    plt::legend();
    string figure1Title;
    figure1Title = "$\\theta$ and $\\dot{\\theta}$ against time";
    plt::title(figure1Title);
    plt::save("Simple pendulum: theta and theta dot against t.svg");
    plt::figure(2);
    plt::plot(theta, thetaDot);
    plt::xlabel("$\\theta$");
    plt::ylabel("$\\dot{\\theta}$");
    plt::title("Phase plot");
    plt::save("Simple pendulum: phase plot of thetaDot against theta.svg");
    return 1;
}
