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
vector<double> SEIR(vector<double> params, double t, vector<double> vars, 
double dt) {
    // Extract parameters
    double a = params[0];
    double beta = params[1];
    double gamma = params[2];
    double delta = params[3];
    double lambda = params[4];
    double mu = params[5];
    double N = params[6];

    // Extract variables
    double S = vars[0];
    double E = vars[1];
    double I = vars[2];
    double R = vars[3];

    // Derivatives and return differentials
    double dSdt = lambda*N - mu*S - beta*I*(1-delta)*S/N;
    double dEdt = beta*I*(1-delta)*S/N - (mu + a) * E;
    double dIdt = a * E - (gamma + mu) * I;
    double dRdt = gamma*I - mu * R;
    return {dSdt*dt, dEdt*dt, dIdt*dt, dRdt*dt};
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
    double S0 = 89.0;
    double E0 = 0.0;
    double I0 = 11.0;
    double R0 = 0.0;
    vector<double> conds = {S0, E0, I0, R0};
    double t0 = 0;
    double tf = 98;

    // Problem parameters
    double a = 1.0;
    double beta = 1.5;
    double gamma = 0.25;
    double delta = 0.0;
    double lambda = 1e-4;
    double mu = 1e-5;
    double N = S0 + E0 + I0 + R0;
    vector<double> params = {a, beta, gamma, delta, lambda, mu, N};

    // Solve problem
    solClass solution = RKF45(SEIR, dtInitial, epsilon, params, t0, tf, conds);
    vector<double> t = solution.t;
    int k = t.size();
    vector<vector<double>> vars = solution.vars;
    vector<double> S = vars[0];
    vector<double> E = vars[1];
    vector<double> I = vars[2];
    vector<double> R = vars[3];

    // Write to file
    ofstream myfile;
    myfile.open("SEIR.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(k), ' ');
    myfile << "t" << std::string(19, ' ');
    myfile << "S" << std::string(19, ' ');
    myfile << "E" << std::string(19, ' ');
    myfile << "I" << std::string(19, ' ');
    myfile << "R" << "\n";
    // Contents
    for (int i = 0 ; i < k; i++) {
        if (i == 0) {
            myfile << i << std::string(1 + (int)log10(k), ' ');
        } else {
            myfile << i << std::string(1 + (int)log10(k) - (int)log10(i), ' ');
        }
        myfile << setprecision(15) << t[i] << " ";
        myfile << setprecision(15) << S[i] << " ";
        myfile << setprecision(15) << E[i] << " ";
        myfile << setprecision(15) << I[i] << " ";
        myfile << setprecision(15) << R[i] << "\n";
    }

    // Plot using matplotlibcpp
    // You will get linting errors for plt::plot, but no build errors if your
    // matplotlibcpp package is installed and set up properly
    plt::figure(1);
    plt::plot(t, S, {{"label", "S"}});
    plt::plot(t, E, {{"label", "E"}});
    plt::plot(t, I, {{"label", "I"}});
    plt::plot(t, R, {{"label", "R"}});
    plt::xlabel("Time (days)");
    plt::legend();
    plt::save("SEIR against time plot.svg");
    return 1;
}
