#include <cmath>
#include <cstdio>

int main() {
    int N = 10;
    double pi = 0;
    for(float k = 0.0; k <= N; k = k + 1) {
         pi = pi + 1/pow(16.0,k)*(4.0/(8.0*k+1.0)-2.0/(8.0*k+4.0)-1.0/(8.0*k+5.0)-1.0/(8.0*k+6.0));
    }

    printf("%.15f", pi);
}
