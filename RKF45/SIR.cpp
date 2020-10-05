#include <RKF45.cpp>

/**
 * Find dt times the RHS of the ODE expressed as a system of first-order equations.
 *
 * @param params   A vector of parameters.
 * @param t        Time value.
 * @param vars     A vector of dependent variable values.
 * @param dt       Step size.
 * @return         vector of dtheta, dthetaDot
 */
std::vector<double> SIR(std::vector<double> params, double t, std::vector<double> vars, double dt) {
    double beta = params[0];
    double gamma = params[1];
    double delta = params[2];
    double N = params[3];

    double S = vars[0];
    double I = vars[1];
    double R = vars[2];

    double dSdt = -beta*I*(1-delta)*S/N;
    double dIdt = beta*I*(1-delta)*S/N - gamma*I;
    double dRdt = gamma*I;

    return {dSdt*dt, dIdt*dt, dRdt*dt};
}

/**
 * @brief          Solves the problem and provides desired output, such as saved plots and data in a textfile.
 */
int main() {
    // Initialize relevant variables
    double epsilon = 1e-11;
    double S0 = 89.0;
    double I0 = 11.0;
    double R0 = 0.0;
    double beta = 1.5;
    double gamma = 0.25;
    double delta = 0.5;
    double N = S0 + I0 + R0;
    std::vector<double> params = {beta, gamma, delta, N};
    double dtInitial = 0.1;
    double t0 = 0;
    double tf = 98;
    // Solve problem
    solClass solution = RKF45(SIR, dtInitial, epsilon, params, t0, tf, {S0, I0, R0});
    std::vector<double> t = solution.t;
    int k = t.size();
    std::vector<std::vector<double>> vars = solution.vars;
    std::vector<double> S = vars[0];
    std::vector<double> I = vars[1];
    std::vector<double> R = vars[2];

    // Write to file
    ofstream myfile;
    myfile.open("SIR.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(k), ' ');
    myfile << "t" << std::string(19, ' ');
    myfile << "S" << std::string(19, ' ');
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
        myfile << setprecision(15) << I[i] << " ";
        myfile << setprecision(15) << R[i] << "\n";
    }

    // Plot using matplotlibcpp
    // You will get linting errors for plt::plot, but no build errors if your
    // matplotlibcpp package is installed and set up properly
    plt::figure(1);
    plt::plot(t, S, {{"label", "S"}});
    plt::plot(t, I, {{"label", "I"}});
    plt::plot(t, R, {{"label", "R"}});
    plt::legend();
    plt::save("SIR against time plot.svg");
    return 1;
}