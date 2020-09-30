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
 * @param x        x in the Lorenz system.
 * @param y        y in the Lorenz system.
 * @param z        z in the Lorenz system.
 * @param dt       Step size.
 * @return         vector of dx, dy, dz
 */
std::vector<double> Lorenz(double sigma, double rho, double beta, double t, double x, double y, double z, double dt) {
    double dx = dt*sigma*(y-x);
    double dy = dt*(x*(rho-z)-y);
    double dz = dt*(x*y - beta*z);
    return {dx, dy, dz};
}

/**
 * @brief Create an object with which we can store solution values, namely i, t, x, y and z.
 */
class solClass {
    public:
        int i; /**< The number of elements each contained vector */
        std::vector<double> t; /**< Array of t values. */
        std::vector<double> x; /**< Array of x values. */
        std::vector<double> y; /**< Array of y values. */
        std::vector<double> z; /**< Array of z values. */
};

/**
 * Solve the problem using the Runge-Kutta-Fehlberg method of 4/5th order.
 *
 * @param f             Function that calculates the RHS of the problem expressed as a system.
 * @param dtInitial     Initial step size.
 * @param epsilon       Error tolerance.
 * @param sigma         Problem parameter.
 * @param rho           Problem parameter.
 * @param beta          Problem parameter.
 * @param t0            Starting time of the simulation.
 * @param tf            End time of the simulation.
 * @param x0            Initial x value.
 * @param y0            Initial y value.
 * @param z0            Initial z value.
 * @return              Solution object containing i, t, x, y and z.
 */
solClass RKF45(std::vector<double>(*f)(double, double, double, double, double, double, double, double), double dtInitial, double epsilon, double sigma, double rho, double beta, double t0, double tf, double x0, double y0, double z0) {
    // Initialize vectors
    std::vector<double> t;
    t.push_back(t0);
    std::vector<double> x;
    x.push_back(x0);
    std::vector<double> y;
    y.push_back(y0);
    std::vector<double> z;
    z.push_back(z0);

    // Initialize other variables
    int i = 0;
    double dt = dtInitial;
    solClass solution;

    // Loop over t[i], increasing it by df at each step, until t[i] = tf
    while ( t[i] < tf ) {
        dt = std::min(dt, tf-t[i]);

        // Approximators for the change in x and y
        std::vector<double> K1 = f(sigma, rho, beta, t[i], x[i], y[i], z[i], dt);
        double k1 = K1[0];
        double l1 = K1[1];
        double m1 = K1[2];
        std::vector<double> K2 = f(sigma, rho, beta, t[i] + dt/4, x[i] + k1/4, y[i] + l1/4.0, z[i] + m1/4.0, dt);
        double k2 = K2[0];
        double l2 = K2[1];
        double m2 = K2[2];
        std::vector<double> K3 = f(sigma, rho, beta, t[i] + 3.0*dt/8.0, x[i] + 3.0*k1/32.0 + 9.0*k2/32.0, y[i] + 3.0*l1/32.0 + 9.0*l2/32.0, z[i] + 3.0*m1/32.0 + 9.0*m2/32.0, dt);
        double k3 = K3[0];
        double l3 = K3[1];
        double m3 = K3[2];
        std::vector<double> K4 = f(sigma, rho, beta, t[i]+12.0*dt/13.0, x[i]+1932.0*k1/2197.0-7200.0*k2/2197.0+7296.0*k3/2197.0, y[i]+1932.0*l1/2197.0-7200.0*l2/2197.0+7296.0*l3/2197.0, z[i]+1932.0*m1/2197.0-7200.0*m2/2197.0+7296.0*m3/2197.0, dt);
        double k4 = K4[0];
        double l4 = K4[1];
        double m4 = K4[2];
        std::vector<double> K5 = f(sigma, rho, beta, t[i]+dt, x[i]+439.0*k1/216.0-8.0*k2+3680.0*k3/513.0-845.0*k4/4104.0, y[i]+439.0*l1/216.0-8.0*l2+3680.0*l3/513.0-845.0*l4/4104.0, z[i]+439.0*m1/216.0-8.0*m2+3680.0*m3/513.0-845.0*m4/4104.0, dt);
        double k5 = K5[0];
        double l5 = K5[1];
        double m5 = K5[2];
        std::vector<double> K6 = f(sigma, rho, beta, t[i]+dt/2.0, x[i]-8.0*k1/27.0+2.0*k2-3544.0*k3/2565.0+1859.0*k4/4104.0-11.0*k5/40.0, y[i]-8.0*l1/27.0+2.0*l2-3544.0*l3/2565.0+1859.0*l4/4104.0-11.0*l5/40.0, z[i]-8.0*m1/27.0+2.0*m2-3544.0*m3/2565.0+1859.0*m4/4104.0-11.0*m5/40.0, dt);
        double k6 = K6[0];
        double l6 = K6[1];
        double m6 = K6[2];

        // x1 and y1 are 4th order approx
        // x2 and y2 are 5th order approx
        double x1 = x[i] + 25.0*k1/216.0+1408.0*k3/2565.0+2197.0*k4/4104.0-k5/5.0;
        double y1 = y[i] + 25.0*l1/216.0+1408.0*l3/2565.0+2197.0*l4/4104.0-l5/5.0;
        double z1 = z[i] + 25.0*m1/216.0+1408.0*m3/2565.0+2197.0*m4/4104.0-m5/5.0;
        double x2 = x[i] + 16.0*k1/135.0+6656.0*k3/12825.0+28561.0*k4/56430.0-9.0*k5/50.0+2.0*k6/55.0;
        double y2 = y[i] + 16.0*l1/135.0+6656.0*l3/12825.0+28561.0*l4/56430.0-9.0*l5/50.0+2.0*l6/55.0;
        double z2 = z[i] + 16.0*m1/135.0+6656.0*m3/12825.0+28561.0*m4/56430.0-9.0*m5/50.0+2.0*m6/55.0;

        // Approximate error in x and y and, if necessary, adjust
        // dt to fix the error
        double Rx = abs(x1-x2)/dt;
        double Ry = abs(y1-y2)/dt;
        double Rz = abs(z1-z2)/dt;
        double sx = pow(epsilon/(2.0*Rx), 0.25);
        double sy = pow(epsilon/(2.0*Ry), 0.25);
        double sz = pow(epsilon/(2.0*Rz), 0.25);
        double R = std::max({Rx, Ry, Rz});
        double s = std::min({sx, sy, sz});
        if (R <= epsilon) {
            t.push_back(t[i]+dt);
            x.push_back(x1);
            y.push_back(y1);
            z.push_back(z1);
            dt *= s;
            i++;
        } else {
            dt *= s;
        }
    }

    // Write solution values to solClass object
    solution.i = i;
    solution.t = t;
    solution.x = x;
    solution.y = y;
    solution.z = z;

    // Return a vector of a vector
    return solution;
}

/**
 * @brief          Solves the problem and provides desired output, such as saved plots and data in a textfile.
 */
int main() {
    // Initialize relevant variables
    double epsilon = 1e-9;
    double x0 = 0.1;
    double y0 = 0.1;
    double z0 = 0.1;
    double sigma = 10;
    double rho = 28;
    double beta = 8.0/3.0;
    double dtInitial = 0.1;
    double t0 = 0;
    double tf = 100;

    // Solve problem
    solClass solution = RKF45(Lorenz, dtInitial, epsilon, sigma, rho, beta, t0, tf, x0, y0, z0);
    int k = solution.i;
    std::vector<double> t = solution.t;
    std::vector<double> x = solution.x;
    std::vector<double> y = solution.y;
    std::vector<double> z = solution.z;

    // Write to file
    ofstream myfile;
    myfile.open("RKF45.txt");
    // Headings
    myfile << "i" << std::string(1 + (int)log10(k), ' ');
    myfile << "t" << std::string(19, ' ');
    myfile << "x" << std::string(17, ' ');
    myfile << "y" << std::string(17, ' ');
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
        myfile << setprecision(15) << y[i] << " ";
        myfile << setprecision(15) << z[i] << "\n";
    }

    // Plot using matplotlibcpp
    // You will get linting errors for plt::plot, but no build errors if your
    // matplotlibcpp package is installed and set up properly
    //plt::figure(1);
    plt::plot3(x, y, z);
    // plt::save("Lorenz 3D phase plot.svg");
    plt::show();
    return 0;
}