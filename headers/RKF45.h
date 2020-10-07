#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

/**
 * @brief Create an object with which we can store solution values, namely i, 
 * t, theta and thetaDot.
 */
class solClass {
    public:
        int i; /**< The number of elements contained within the vectors 
        below */
        vector<double> t; /**< Vector of t values. */
        vector<vector<double>> vars; /**< Vector of vectors of dependent 
        variable values. */
};

/**
 * Multiply a vector by a scalar
 * 
 * @param vec      Vector to be multiplied.
 * @param scalar   A floating point number the vector is to be multiplied by.
 * @return         vec * scalar
 */
vector<double> vecMult(vector<double> vec, double scalar) {
    int N = vec.size();
    vector<double> multvec = {};
    for (int i = 0; i < N; i++) {
        multvec.push_back(vec[i]*scalar);
    }
    return multvec;
}

/**
 * Add vectors within the vector of vectors vvec
 * 
 * @param vvec          Vector of vectors that are to be added.
 * @return              Vector that is the sum of the vectors contained within
 * vvec.
 */
vector<double> vecAdd(vector<vector<double>> vvec) {
    int N = vvec.size();
    int M = (vvec[0]).size();
    vector<double> addedVec(M, 0);
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++) {
            addedVec[j] += vvec[i][j];
        }
    }
    return addedVec;
}

/**
 * Find the absolute value of a specified vector
 * 
 * @param vec           Vector whose absolute value is to be determined.
 * @return              abs(vec)
 */
vector<double> vecAbs(vector<double> vec) {
    int N = vec.size();
    vector<double> vecAbs;
    for (int i = 0; i < N ; i++) {
        vecAbs.push_back(abs(vec[i]));
    }
    return vecAbs;
}

/**
 * Find the maximum value of a specified vector
 *
 * @param vec           Vector whose maximum value is being found.
 * @return              Maximum element of vec.
 */
double vecMax(vector<double> vec) {
    double max = *std::max_element(vec.begin(), vec.end());
    return max;
}

/**
 * Solve the problem using the Runge-Kutta-Fehlberg method of 4/5th order.
 *
 * @param f             Function that calculates the RHS of the problem 
 * expressed as a system.
 * @param dtInitial     Initial step size.
 * @param epsilon       Error tolerance.
 * @param g             Acceleration due to gravity in metres per second 
 * squared.
 * @param l             Length of the pendulum rod in metres.
 * @param t0            Start time of the simulation.
 * @param tf            End time of the simulation.
 * @param theta0        Initial angle from the positive x-axis in radians.
 * @param thetaDot0     Rate of change of theta with respect to time at t0.
 * @return              Solution object containing i, t, theta and thetaDot.
 */
solClass RKF45(vector<double>(*f)(vector<double>, double, 
vector<double>, double), double dtInitial, double epsilon, 
vector<double> params, double t0, double tf, vector<double> conds) {
    // Initialize vectors
    vector<double> t;
    t.push_back(t0);
    vector<vector<double>> vars;
    vars.push_back(conds);

    // Initialize other variables
    int i = 0;
    double dt = dtInitial;
    solClass solution;

    // Loop over t[i], increasing it by df at each step, until t[i] = tf
    while (t[i] < tf ) {
        dt = std::min(dt, tf-t[i]);

        // Approximators for the change in theta and thetaDot
        // Blank lines in between to help with readability;
        // Split over multiple lines to prevent lines exceeding 80 chars
        vector<double> K1 = f(params, t[i], vars[i], dt);

        vector<double> K2 = f(params, t[i] + dt/4, vecAdd({vars[i], 
        vecMult(K1, 0.25)}), dt);

        vector<double> K3 = f(params, t[i] + 3.0*dt/8.0, vecAdd({vars[i], 
        vecMult(K1, 3.0/32.0), vecMult(K2, 9.0/32.0)}), dt);

        vector<double> K4 = f(params, t[i]+12.0*dt/13.0, vecAdd({vars[i], 
        vecMult(K1, 1932.0/2197.0), vecMult(K2, -7200.0/2197.0), 
        vecMult(K3, 7296.0/2197.0)}), dt);

        vector<double> K5 = f(params, t[i]+dt, vecAdd({vars[i], 
        vecMult(K1, 439.0/216.0), vecMult(K2, -8.0), 
        vecMult(K3, 3680.0/513.0), vecMult(K4, -845.0/4104.0)}), dt);

        vector<double> K6 = f(params, t[i]+dt/2.0, 
        vecAdd({vars[i], vecMult(K1, -8.0/27.0), vecMult(K2, 2.0), 
        vecMult(K3, -3544.0/2565.0), vecMult(K4, 1859.0/4104.0), 
        vecMult(K5, -11.0/40.0)}), dt);

        // theta1 and thetaDot1 are 4th order approx
        // theta2 and thetaDot2 are 5th order approx
        vector<double> vars1 = vecAdd({vars[i], vecMult(K1, 25.0/216.0), 
        vecMult(K3, 1408.0/2565.0), vecMult(K4, 2197.0/4104.0), 
        vecMult(K5, -1.0/5.0)});

        vector<double> vars2 = vecAdd({vars[i], vecMult(K1, 16.0/135.0), 
        vecMult(K3, 6656.0/12825.0), vecMult(K4, 28561.0/56430.0), 
        vecMult(K5, -9.0/50.0), vecMult(K6, 2.0/55.0)});
        // Approximate error in theta and thetaDot and, if necessary, adjust
        // dt to fix the error
        vector<double> Rvars = vecMult(vecAbs(vecAdd(
            {vars1, vecMult(vars2, -1.0)})), 1.0/dt);
        double R = vecMax(Rvars);
        double s = pow((epsilon/(2.0*R)), 0.25);
        if (R <= epsilon) {
            t.push_back(t[i]+dt);
            vars.push_back(vars1);
            i++;
        }
        dt *= s;
    }

    // Transpose vars to make it easier to obtain solution values
    int N1 = vars.size();
    int N2 = vars[0].size();
    vector<vector<double>> varsTrans(N2, vector<double>(N1));
    for (int l = 0; l < N1; l++) {
        for (int m = 0; m < N2 ; m++) {
            varsTrans[m][l] = vars[l][m];
        }
    }

    // Write solution values to solClass object
    solution.i = i;
    solution.t = t;
    solution.vars = varsTrans;

    // Return a vector of a vector
    return solution;
}