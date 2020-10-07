#include <integrate.h>
#include <iostream>

double f(double x, std::vector<double> params) {
    return pow(x, 2);
}

int main() {
    double integralAd = adaptive_simpsons(0, 10, 1e-10, {0.0}, f);
    double intSimp = simpsons(0, 10, 1e3, {0.0}, f);
    // std::cout << integralAd;
    std::cout << intSimp;
}