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

    // Calculate period as double the integral from thetaMin to thetaMax
    return integral;
}