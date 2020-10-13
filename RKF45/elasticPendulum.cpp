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
 * @return         vector of differentials
 */
vector<double> doubPen(vector<double> params, double t, vector<double> vars, 
double dt) {
    // Extract parameters
    double g = params[0];
    double l0 = params[1];
    double k = params[2];
    double m = params[3];

    // Extract variables from vars
    double x = vars[0];
    double xDot = vars[1];
    double theta = vars[2];
    double thetaDot = vars[3];

    // Derivatives
    double xDDot = (l0+x)*pow(thetaDot, 2) - k*x/m + g * cos(theta);
    double thetaDDot = (-g/(l0+x))*sin(theta)-(2*xDot/(l0+x))*thetaDot;

    // Differentials
    double dx = dt * xDot;
    double dxDot = dt * xDDot;
    double dTheta = dt * thetaDot;
    double dThetaDot = dt * thetaDDot;
    return {dx, dxDot, dTheta, dThetaDot};
}

/**
 * @brief Solves the problem and provides desired output, such as saved plots
 * and data in a textfile.
 */
int main() {
    // Initial conditions and domain of integration
    double x0 = 1.0;
    double xDot0 = 0.0;
    double theta0 = M_PI/2;
    double thetaDot0 = 0.0;
    vector<double> conds = {x0, xDot0, theta0, thetaDot0};
    double t0 = 0.0;
    double tf = 20.0;

    // Problem parameters
    double g = 9.81;
    double l0 = 1.0;
    double k = 1.0;
    double m = 1.0;
    vector<double> params = {g, l0, k, m};

    // Other parameters
    double epsilon = 1e-9;
    double dtInitial = 0.1;

    // Solve problem
    solClass solution = RKF45(doubPen, dtInitial, epsilon, params, t0, tf, 
    conds);
    vector<double> t = solution.t;
    vector<vector<double>> vars = solution.vars;
    int K = t.size();

    // Extract solution values from vars
    vector<double> x = vars[0];
    vector<double> xDot = vars[1];
    vector<double> theta = vars[2];
    vector<double> thetaDot = vars[3];

    // Write to file
    ofstream myfile;
    myfile.open("elasticPendulum.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(K), ' ');
    myfile << "t" << std::string(19, ' ');
    myfile << "x" << std::string(17, ' ');
    myfile << "xDot" << std::string(16, ' ');
    myfile << "theta" << std::string(13, ' ');
    myfile << "thetaDot" << "\n";
    // Contents
    for (int i = 0 ; i < K; i++) {
        // Spacing in file between columns
        if (i == 0) {
            myfile << i << std::string(1 + (int)log10(K), ' ');
        } else {
            myfile << i << std::string(1 + (int)log10(K) - (int)log10(i), ' ');
        }

        // Table entries in file
        myfile << setprecision(15) << t[i] << " ";
        myfile << setprecision(15) << x[i] << " ";
        myfile << setprecision(15) << xDot[i] << " ";
        myfile << setprecision(15) << theta[i] << " ";
        myfile << setprecision(15) << thetaDot[i] << "\n";
    }

    // Plot using matplotlibcpp
    // You will get linting errors for plt::plot, but no build errors if your
    // matplotlibcpp package is installed and set up properly
    plt::figure(1);
    plt::plot(t, x, {{"label", "$x$"}});
    plt::plot(t, xDot, {{"label", "$\\dot{x}$"}});
    plt::xlabel("$t$");
    plt::legend();
    string figure1Title;
    figure1Title = "$x$ and $\\dot{x}$ against time";
    plt::title(figure1Title);
    plt::save("Elastic pendulum: x and xDot against time plot.svg");
    plt::figure(2);
    plt::plot(t, theta, {{"label", "$\\theta$"}});
    plt::plot(t, thetaDot, {{"label", "$\\dot{\\theta}$"}});
    plt::xlabel("$t$");
    plt::legend();
    string figure2Title;
    figure2Title = "$\\theta$ and $\\dot{\\theta}$ against time";
    plt::title(figure2Title);
    plt::save("Elastic pendulum: theta and thetaDot against time plot.svg");
    plt::figure(3);
    plt::plot(theta, thetaDot);
    plt::xlabel("$\\theta$");
    plt::ylabel("$\\dot{\\theta}$");
    plt::title("Phase plot for theta");
    plt::save("Elastic pendulum: phase plot of theta dot against theta.svg");
    plt::figure(4);
    plt::plot(x, xDot);
    plt::xlabel("$x$");
    plt::ylabel("$\\dot{x}$");
    plt::title("Phase plot for x");
    plt::save("Elastic pendulum: phase plot of ptheta2 against theta2.svg");
    return 1;
}
