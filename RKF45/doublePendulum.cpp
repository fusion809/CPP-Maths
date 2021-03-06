// From matplotlib-cpp-git AUR package
#include <matplotlib-cpp/matplotlibcpp.h>
namespace plt = matplotlibcpp;
#include <iostream>
// Required for setprecision
#include <iomanip>
#include <string>
#include <fstream>
// RKF45
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
vector<double> doubPen(vector<double> params, double t, 
vector<double> vars, double dt) {
    // Extract parameters
    double g  = params[0];
    double l1 = params[1];
    double l2 = params[2];
    double m1 = params[3];
    double m2 = params[4];

    // Extract variables from vars
    double theta1  = vars[0];
    double pTheta1 = vars[1];
    double theta2  = vars[2];
    double pTheta2 = vars[3];

    // Define variables used to simplify the ODE system
    double C1 = pTheta1*pTheta2*sin(theta1-theta2)/(l1*l2*(m1+m2*
    (pow(sin(theta1-theta2), 2))));
    double C2 = ((pow(l2,2))*m2*pow(pTheta1, 2) + (pow(l1, 2))*(m1+m2)*
    pow(pTheta2, 2) - l1*l2*m2*pTheta1*pTheta2*cos(theta1-theta2))/
    (2*(pow(l1,2))*(pow(l2,2))*(m1+m2*pow(sin(theta1-theta2), 2)))*
    sin(2*(theta1-theta2));

    // Derivatives
    double thetaDot1 = (l2 * pTheta1 - l1*pTheta2 * cos(theta1))/((pow(l1, 2))*
    l2*(m1+m2*((pow(sin(theta1-theta2), 2)))));
    double pThetaDot1 = -(m1+m2)*g*l1*sin(theta1) - C1 + C2;
    double thetaDot2 = (l1*(m1+m2)*pTheta2-l2*m2*pTheta1*cos(theta1-theta2))/
    (l1*(pow(l2, 2))*m2*(m1+m2*(pow(sin(theta1-theta2), 2))));
    double pThetaDot2 = -m2*g*l2*sin(theta2)+C1-C2;

    // Differentials and return them
    double dTheta1 = dt * thetaDot1;
    double dPTheta1 = dt * pThetaDot1;
    double dTheta2 = dt * thetaDot2;
    double dPTheta2 = dt * pThetaDot2;
    return {dTheta1, dPTheta1, dTheta2, dPTheta2};
}

/**
 * @brief Solves the problem and provides desired output, such as saved plots 
 * and data in a textfile.
 */
int main() {
    // Initial conditions and domain of integration
    double theta10 = M_PI/2;
    double pTheta10 = 0;
    double theta20 = M_PI/2;
    double pTheta20 = 0;
    vector<double> conds = {theta10, pTheta10, theta20, pTheta20};
    double t0 = 0.0;
    double tf = 100.0;

    // Problem parameters
    double g = 9.81;
    double l1 = 1.0;
    double l2 = 1.0;
    double m1 = 1.0;
    double m2 = 1.0;
    vector<double> params = {g, l1, l2, m1, m2};

    // Other parameters
    double epsilon = 1e-9;
    double dtInitial = 0.1;

    // Solve problem
    auto start = high_resolution_clock::now(); 
    solClass solution = RKF45(doubPen, dtInitial, epsilon, params, t0, tf,
    conds);
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<milliseconds>(stop - start);
    cout << "RKF45 took " << duration.count() << " ms to run" << endl;
    vector<double> t = solution.t;
    vector<vector<double>> vars = solution.vars;
    int k = t.size();

    // Extract solution values from vars
    vector<double> theta1 = vars[0];
    vector<double> pTheta1 = vars[1];
    vector<double> theta2 = vars[2];
    vector<double> pTheta2 = vars[3];

    // Initialize pendulum coordinates
    vector<double> x1(k, 0.0);
    vector<double> x2(k, 0.0);
    vector<double> y1(k, 0.0);
    vector<double> y2(k, 0.0);

    // Write to file
    ofstream myfile;
    myfile.open("doublePendulum.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(k), ' ');
    myfile << "t" << std::string(19, ' ');
    myfile << "theta1" << std::string(12, ' ');
    myfile << "pTheta1" << std::string(12, ' ');
    myfile << "theta2" << std::string(12, ' ');
    myfile << "pTheta2" << "\n";
    // Contents
    for (int i = 0 ; i < k; i++) {
        // Spacing in file between columns
        if (i == 0) {
            myfile << i << std::string(1 + (int)log10(k), ' ');
        } else {
            myfile << i << std::string(1 + (int)log10(k) - (int)log10(i), ' ');
        }

        // Calculate pendulum Cartesian coordinates
        x1[i] = l1*sin(theta1[i]);
        y1[i] = -l1*cos(theta1[i]);
        x2[i] = x1[i] + l2*sin(theta2[i]);
        y2[i] = y1[i] - l2*cos(theta2[i]);

        // Table entries in file
        myfile << setprecision(15) << t[i] << " ";
        myfile << setprecision(15) << theta1[i] << " ";
        myfile << setprecision(15) << pTheta1[i] << " ";
        myfile << setprecision(15) << theta2[i] << " ";
        myfile << setprecision(15) << pTheta2[i] << "\n";
    }

    // Plot using matplotlibcpp
    // You will get linting errors for plt::plot, but no build errors if your
    // matplotlibcpp package is installed and set up properly
    plt::figure(1);
    plt::plot(t, theta1, {{"label", "$\\theta_1$"}});
    plt::plot(t, pTheta1, {{"label", "$p_{\\theta_1}$"}});
    plt::plot(t, theta2, {{"label", "$\\theta_2$"}});
    plt::plot(t, pTheta2, {{"label", "$p_{\\theta_2}$"}});
    plt::xlabel("$t$");
    plt::legend();
    string fig1Title;
    fig1Title = "$\\theta_1$, $p_{\\theta_1}$, $\\theta_2$, and ";
    fig1Title += "$p_{\\theta_2}$ against time";
    plt::title(fig1Title);
    string fig1Filename = "Double pendulum: theta1, pTheta1, theta2 and";
    fig1Filename += "pTheta2 against time plot.svg";
    plt::save(fig1Filename);
    plt::figure(2);
    plt::plot(theta1, pTheta1);
    plt::xlabel("$\\theta_1$");
    plt::ylabel("$p_{\\theta_1}$");
    plt::title("Phase plot");
    plt::save("Double pendulum: phase plot of pTheta1 against theta1.svg");
    plt::figure(3);
    plt::plot(theta2, pTheta2);
    plt::xlabel("$\\theta_2$");
    plt::ylabel("$p_{\\theta_2}$");
    plt::title("Phase plot");
    plt::save("Double pendulum: phase plot of pTheta2 against theta2.svg");
    plt::figure(4);
    plt::plot(theta1, theta2);
    plt::xlabel("$\\theta_1$");
    plt::ylabel("$\\theta_2$");
    plt::title("Phase plot");
    plt::save("Double pendulum: phase plot of theta2 against theta1.svg");
    plt::figure(5);
    plt::plot(x1, y1, {{"label", "Pendulum 1 path"}});
    plt::plot(x2, y2, {{"label", "Pendulum 2 path"}});
    plt::legend();
    plt::title("Path of the double pendulum");
    plt::save("Double pendulum path.svg");
    return 1;
}
