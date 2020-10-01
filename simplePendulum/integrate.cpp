#include <integrate.h>
#include <iostream>
#include <chrono>
using namespace std::chrono;
using namespace std;

/**
 * Calculate the inverse of the time derivative of theta.
 * 
 * @param theta    Angle from the positive x-axis.
 * @param params   A vector containing all parameters of the problem, namely:
 * g, l, theta0 and thetaDot0.
 * @return         1/thetaDot.
 */
double invThetaDot(double theta, std::vector<double> params) {
    // Extract parameters. 
    double g = params[0];
    double l = params[1];
    double theta0 = params[2];
    double thetaDot0 = params[3];

    return 1/sqrt(pow(thetaDot0, 2.0) + 2*g/l * (sin(theta0)-sin(theta)));
}

int main() {
    // Start timer
    auto start = high_resolution_clock::now(); 
    double intSimps = simpsons(-1e-6, -M_PI + 1e-6, 1e8, {9.81, 1.0, 0.0, 0.0}, invThetaDot);
    // End timer
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<milliseconds>(end - start); 
    cout << "intSimps is " << intSimps << ". Duration of the calculation is " << duration.count()/1e3 << " seconds" << endl;
}

// Julia took 7.48s to integrate 1/sqrt(-19.62*sin(x)) from -1e-5 to -pi+1e-5 with N=1e8, whereas C++ took 18.18s.