#include "RKF45.h"

/**
 * Find dt times the RHS of the ODE expressed as a system of first-order equations.
 *
 * @param params   A vector of parameters.
 * @param t        Time value.
 * @param vars     A vector of dependent variable values.
 * @param dt       Step size.
 * @return         vector of differentials
 */
std::vector<double> Chen(std::vector<double> params, double t, std::vector<double> vars, double dt) {
    // Extract parameter values
    double a = params[0];
    double b = params[1];
    double c = params[2];

    // Extract variable values
    double x = vars[0];
    double y = vars[1];
    double z = vars[2];

    // Derivatives and differentials
    double dxdt = a*(y-x);
    double dydt = (c-a)*x-x*z+c*y;
    double dzdt = x*y-b*z;
    return {dt*dxdt, dt*dydt, dt*dzdt};
}

/**
 * @brief          Solves the problem and provides desired output, such as saved plots and data in a textfile.
 */
int main() {
    // Solution parameters
    double dtInitial = 0.1;
    double epsilon = 1e-9;

    // Initial conditions and domain of integration
    double x0 = -0.1;
    double y0 = 0.5;
    double z0 = -0.6;
    std::vector<double> conds = {x0, y0, z0};
    double t0 = 0;
    double tf = 60;

    // Problem parameters
    double a = 40;
    double b = 3;
    double c = 28;
    std::vector<double> params = {a, b, c};

    // Solve problem
    solClass solution = RKF45(Chen, dtInitial, epsilon, params, t0, tf, conds);
    std::vector<double> t = solution.t;
    int k = t.size();
    std::vector<std::vector<double>> vars = solution.vars;
    std::vector<double> x = vars[0];
    std::vector<double> y = vars[1];
    std::vector<double> z = vars[2];

    // Write to file
    ofstream myfile;
    myfile.open("Chen.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(k), ' ');
    myfile << "t" << std::string(19, ' ');
    myfile << "x" << std::string(19, ' ');
    myfile << "y" << std::string(19, ' ');
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
    plt::save("Chen 3D phaseplot.svg");
    plt::figure(3);
    plt::clf();
    plt::plot(t, x, {{"label", "$x$"}});
    plt::plot(t, y, {{"label", "$y$"}});
    plt::plot(t, z, {{"label", "$z$"}});
    plt::legend();
    plt::save("Chen x, y and z against t plot.svg");
    return 1;
}
