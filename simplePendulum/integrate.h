#include <cmath>
#include <vector>

// Use Chebyshev-Gauss quadrature
/**
 * @brief Use Chebyshev-Gauss quadrature to integrate f from a to b with N nodes.
 * 
 * @param a        Start of the domain of integration.
 * @param b        End of the domain of integration.
 * @param N        Number of nodes used.
 * @param params   A vector containing all parameters required by f.
 * @param f        The function be integrated.
 * @return         Approximated value of the integral of f from a to b.
 */
double chebGaussQuad(double a, double b, double N, std::vector<double> params, double(*f)(double, std::vector<double> params)) {
    // Initialize variables
    double nodes = 0;
    double integrand = 0;
    double transformedGrid = 0;
    double integral = 0;

    // Loop over each node
    for (int i = 1; i < N +1 ; i++) {
        nodes = cos((2.0*i-1.0)*M_PI/(2.0*N));
        transformedGrid = (b-a)*nodes/2.0 + (a+b)/2.0;
        integrand = sqrt(1-pow(nodes, 2.0))*f(transformedGrid, params);
        integral += ((b-a)/2.0) * (M_PI/N) * integrand;
    }

    // Return integral
    return integral;
};

/**
 * @brief Uses Simpson's rule with N steps to integrate f(x, param)
 * 
 * @param a        Start of the domain of integration.
 * @param b        End of the domain of integration.
 * @param N        Step number.
 * @param params   A vector containing all parameters for function f.
 * @param f        A funcction of x and params. 
 * @return         The integral of f(x, param) from a to b.
 */
double simpsons(double a, double b, int N, std::vector<double> params, double(*f)(double, std::vector<double>)) {
    // Step size
    if (N % 2 == 1) {
        return 0;
    }
    double step = (b-a)/N;
    double grid = a;
    double integral = 0;

    for (int i = 0; i <= N ; i++) {
        if ( ( i == 0 ) || (i == N) ) {
            integral += (step/3) * f(grid, params);
        } else if ( (i - 1) % 2 == 0) {
            integral += (2*step/3) * f(grid, params);
        } else {
            integral += (4*step/3) * f(grid, params);
        }
        grid += step;
    }

    return integral;
};