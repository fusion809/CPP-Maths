#include <RKF45.cpp>

/**
 * Find dt times the RHS of the ODE expressed as a system of first-order equations.
 *
 * @param g        Acceleration due to gravity in metres per second squared.
 * @param l        Length of the pendulum rod in metres.
 * @param t        Time value. Largely included for the sake of generality, not actually used.
 * @param theta    Angle from the positive x-axis (positive values = above the x-axis).
 * @param thetaDot Rate of change of theta with respect to time.
 * @param dt       Step size.
 * @return         vector of dtheta, dthetaDot
 */
std::vector<double> Lorenz(std::vector<double> params, double t, std::vector<double> vars, double dt) {
    double sigma = params[0];
    double beta = params[1];
    double rho = params[2];
    double x = vars[0];
    double y = vars[1];
    double z = vars[2];
    double dxdt = sigma*(y-x);
    double dydt = x*(rho-z) - y;
    double dzdt = x*y - beta*z;
    return {dt*dxdt, dt*dydt, dt*dzdt};
}

/**
 * @brief          Solves the problem and provides desired output, such as saved plots and data in a textfile.
 */
int main() {
    // Initialize relevant variables
    double epsilon = 1e-9;
    double x0 = 1.0;
    double y0 = 1.0;
    double z0 = 1.0;
    double sigma = 10.0;
    double beta = 8.0/3.0;
    double rho = 28.0;
    std::vector<double> params = {sigma, beta, rho};
    double dtInitial = 0.1;
    double t0 = 0;
    double tf = 60;
    // Solve problem
    solClass solution = RKF45(Lorenz, dtInitial, epsilon, params, t0, tf, {x0, y0, z0});
    std::vector<double> t = solution.t;
    int k = t.size();
    std::vector<std::vector<double>> vars = solution.vars;
    std::vector<double> x = vars[0];
    std::vector<double> y = vars[1];
    std::vector<double> z = vars[2];

    // Write to file
    ofstream myfile;
    myfile.open("Lorenz.txt");
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
    plt::save("Lorenz 3D phaseplot.svg");
    plt::figure(3);
    plt::clf();
    plt::plot(t, x);
    plt::plot(t, y);
    plt::plot(t, z);
    plt::save("x, y, z against t plot.svg");
    return 1;
}
