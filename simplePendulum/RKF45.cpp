#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <utility>
#include <unistd.h>
#include <algorithm>
// Required for setprecision
#include <iomanip>
// From matplotlib-cpp-git AUR package
#include <matplotlib-cpp/matplotlibcpp.h>
using namespace std;
namespace plt = matplotlibcpp;

/**
 * Find the RHS of the ODE expressed as a system of first-order equations.
 *
 * @param g        Acceleration due to gravity in metres per second squared.
 * @param l        Length of the pendulum rod in metres.
 * @param t        Time value. Largely included for the sake of generality, not actually used.
 * @param theta    Angle from the positive x-axis (positive values = above the x-axis).
 * @param thetaDot Rate of change of theta with respect to time.
 * @param dt       Step size.
 * @return         vector of dtheta, dthetaDot
 */
std::vector<double> simpPen(double g, double l, double t, double theta, double thetaDot, double dt) {
    double thetaDDot = -g/l * cos(theta);
    return {dt*thetaDot, dt*thetaDDot};
}

/**
 * @brief Create an object with which we can store solution values, namely i, t, theta and thetaDot.
 */
class solClass {
    public:
        int i; /**< The number of elements contained within the vectors below */
        std::vector<double> t; /**< Vector of t values. */
        std::vector<double> theta; /**< Vector of theta values. */
        std::vector<double> thetaDot; /**< Vector of theta dot values. */
};

/**
 * Calculate the period of the problem.
 *
 * @param g             Acceleration due to gravity in metres per second squared.
 * @param l             Length of the pendulum rod in metres.
 * @param N             Number of grid points used to calculate the period.
 * @param theta0        Initial angle from the positive x-axis.
 * @param thetaDot0     Initial rate of change of theta with respect to t.
 * @return              Period of the problem in seconds.
 */
double periodCalc(double g, double l, int N, double theta0, double thetaDot0) {
    // Initialize relevant variables
    double thetaMin = asin((pow(thetaDot0,2.0) + 2.0*(g/l)*sin(theta0))/(2*g/l));
    double thetaMax = -thetaMin - M_PI;
    double nodes = 0;
    double integrand = 0;
    double transformedGrid = 0;
    double integral = 0;

    // Loop over each node
    for (int i = 1; i < N +1 ; i++) {
        nodes = cos((2.0*i-1.0)*M_PI/(2.0*N));
        transformedGrid = (thetaMax-thetaMin)*nodes/2.0 + (thetaMax + thetaMin)/2.0;
        integrand = sqrt(1-pow(nodes, 2.0))*pow(pow(thetaDot0, 2.0)+2.0*g/l * (sin(theta0)-sin(transformedGrid)), -0.5);
        integral += ((thetaMax-thetaMin)/2.0) * (M_PI/N) * integrand;
    }

    // Calculate period as double the integral from thetaMin to thetaMax
    double period = 2.0*abs(integral);
    return period;
}

/**
 * Solve the problem using the Runge-Kutta-Fehlberg method of 4/5th order.
 *
 * @param f             Function that calculates the RHS of the problem expressed as a system.
 * @param dtInitial     Initial step size.
 * @param epsilon       Error tolerance.
 * @param g             Acceleration due to gravity in metres per second squared.
 * @param l             Length of the pendulum rod in metres.
 * @param t0            Start time of the simulation.
 * @param tf            End time of the simulation.
 * @param theta0        Initial angle from the positive x-axis in radians.
 * @param thetaDot0     Rate of change of theta with respect to time at t0.
 * @return              Solution object containing i, t, theta and thetaDot.
 */
solClass RKF45(std::vector<double>(*f)(double, double, double, double, double, double), double dtInitial, double epsilon, double g, double l, double t0, double tf, double theta0, double thetaDot0) {
    // Initialize vectors
    std::vector<double> t;
    t.push_back(t0);
    std::vector<double> theta;
    theta.push_back(theta0);
    std::vector<double> thetaDot;
    thetaDot.push_back(thetaDot0);

    // Initialize other variables
    int i = 0;
    double dt = dtInitial;
    solClass solution;

    // Loop over t[i], increasing it by df at each step, until t[i] = tf
    while ( t[i] < tf ) {
        dt = std::min(dt, tf-t[i]);

        // Approximators for the change in theta and thetaDot
        std::vector<double> K1 = f(g, l, t[i], theta[i], thetaDot[i], dt);
        double k1 = K1[0];
        double l1 = K1[1];
        std::vector<double> K2 = f(g, l, t[i] + dt/4, theta[i] + k1/4, thetaDot[i] + l1/4.0, dt);
        double k2 = K2[0];
        double l2 = K2[1];
        std::vector<double> K3 = f(g, l, t[i] + 3.0*dt/8.0, theta[i] + 3.0*k1/32.0 + 9.0*k2/32.0, thetaDot[i] + 3.0*l1/32.0 + 9.0*l2/32.0, dt);
        double k3 = K3[0];
        double l3 = K3[1];
        std::vector<double> K4 = f(g, l, t[i]+12.0*dt/13.0, theta[i]+1932.0*k1/2197.0-7200.0*k2/2197.0+7296.0*k3/2197.0, thetaDot[i]+1932.0*l1/2197.0-7200.0*l2/2197.0+7296.0*l3/2197.0, dt);
        double k4 = K4[0];
        double l4 = K4[1];
        std::vector<double> K5 = f(g, l, t[i]+dt, theta[i]+439.0*k1/216.0-8.0*k2+3680.0*k3/513.0-845.0*k4/4104.0, thetaDot[i]+439.0*l1/216.0-8.0*l2+3680.0*l3/513.0-845.0*l4/4104.0, dt);
        double k5 = K5[0];
        double l5 = K5[1];
        std::vector<double> K6 = f(g, l, t[i]+dt/2.0, theta[i]-8.0*k1/27.0+2.0*k2-3544.0*k3/2565.0+1859.0*k4/4104.0-11.0*k5/40.0, thetaDot[i]-8.0*l1/27.0+2.0*l2-3544.0*l3/2565.0+1859.0*l4/4104.0-11.0*l5/40.0, dt);
        double k6 = K6[0];
        double l6 = K6[1];

        // theta1 and thetaDot1 are 4th order approx
        // theta2 and thetaDot2 are 5th order approx
        double theta1 = theta[i] + 25.0*k1/216.0+1408.0*k3/2565.0+2197.0*k4/4104.0-k5/5.0;
        double thetaDot1 = thetaDot[i] + 25.0*l1/216.0+1408.0*l3/2565.0+2197.0*l4/4104.0-l5/5.0;
        double theta2 = theta[i] + 16.0*k1/135.0+6656.0*k3/12825.0+28561.0*k4/56430.0-9.0*k5/50.0+2.0*k6/55.0;
        double thetaDot2 = thetaDot[i] + 16.0*l1/135.0+6656.0*l3/12825.0+28561.0*l4/56430.0-9.0*l5/50.0+2.0*l6/55.0;

        // Approximate error in theta and thetaDot and, if necessary, adjust
        // dt to fix the error
        double RTheta = abs(theta1-theta2)/dt;
        double RThetaDot = abs(thetaDot1-thetaDot2)/dt;
        double sTheta = pow(epsilon/(2.0*RTheta), 0.25);
        double sThetaDot = pow(epsilon/(2.0*RThetaDot), 0.25);
        double R = max(RTheta, RThetaDot);
        double s = min(sTheta, sThetaDot);
        if (R <= epsilon) {
            t.push_back(t[i]+dt);
            theta.push_back(theta1);
            thetaDot.push_back(thetaDot1);
            dt *= s;
            i++;
        } else {
            dt *= s;
        }
    }

    // Write solution values to solClass object
    solution.i = i;
    solution.t = t;
    solution.theta = theta;
    solution.thetaDot = thetaDot;

    // Return a vector of a vector
    return solution;
}

/**
 * @brief          Solves the problem and provides desired output, such as saved plots and data in a textfile.
 */
int main() {
    // Initialize relevant variables
    int N = 1e3;
    double epsilon = 1e-12;
    double theta0 = 0;
    double thetaDot0 = 0;
    double g = 9.81;
    double l = 1.0;
    double dtInitial = 0.1;
    double t0 = 0;
    double tf = t0 + 4*periodCalc(g, l, N, theta0, thetaDot0);

    // Solve problem
    solClass solution = RKF45(simpPen, dtInitial, epsilon, g, l, t0, tf, theta0, thetaDot0);
    int k = solution.i;
    std::vector<double> t = solution.t;
    std::vector<double> theta = solution.theta;
    std::vector<double> thetaDot = solution.thetaDot;

    // Write to file
    ofstream myfile;
    myfile.open("RKF45.txt");
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
    plt::title("theta against time");
    plt::save("theta against t.svg");
    plt::figure(2);
    plt::plot(t, thetaDot);
    plt::xlabel("t");
    plt::ylabel("thetaDot");
    plt::title("thetaDot against time");
    plt::save("thetaDot against t.svg");
    plt::figure(3);
    plt::plot(theta, thetaDot);
    plt::xlabel("theta");
    plt::ylabel("thetaDot");
    plt::title("Phase plot");
    plt::save("Phase plot of thetaDot against theta.svg");
    return 1;
}
