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
std::vector<double> simpPen(std::vector<double> params, double t, std::vector<double> vars, double dt) {
    double g = params[0];
    double l = params[1];
    double theta = vars[0];
    double thetaDot = vars[1];
    double thetaDDot = -g/l * cos(theta);
    return {dt*thetaDot, dt*thetaDDot};
}

/**
 * @brief          Solves the problem and provides desired output, such as saved plots and data in a textfile.
 */
int main() {
    // Initialize relevant variables
    double epsilon = 1e-11;
    double theta0 = 0;
    double thetaDot0 = 0;
    double g = 9.81;
    double l = 1.0;
    std::vector<double> params = {g, l};
    double dtInitial = 0.1;
    double t0 = 0;
    double tf = 10;
    // Solve problem
    solClass solution = RKF45(simpPen, dtInitial, epsilon, params, t0, tf, {theta0, thetaDot0});
    std::vector<double> t = solution.t;
    int k = t.size();
    std::vector<std::vector<double>> vars = solution.vars;
    std::vector<double> theta = vars[0];
    std::vector<double> thetaDot = vars[1];

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
    plt::plot(t, theta);
    plt::xlabel("t");
    plt::ylabel("theta");
    string figure1Title;
    figure1Title = "theta against time";
    plt::title(figure1Title);
    plt::save("theta against t.svg");
    plt::figure(2);
    plt::plot(t, thetaDot);
    plt::xlabel("t");
    plt::ylabel("thetaDot");
    string figure2Title;
    figure2Title = "thetaDot against time";
    plt::title(figure2Title);
    plt::save("thetaDot against t.svg");
    plt::figure(3);
    plt::plot(theta, thetaDot);
    plt::xlabel("theta");
    plt::ylabel("thetaDot");
    plt::title("Phase plot");
    plt::save("Phase plot of thetaDot against theta.svg");
    return 1;
}
