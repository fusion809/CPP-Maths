// This program is meant to solve the Sturm-Liouville equation
#include <algorithm>
#include <cmath>
#include <math.h>
#include <iostream>
#include <numeric>
#include <vector>
#include <fstream>
#undef NDEBUG
// The following include may need to be adjusted on non-Gentoo Linux platforms
#include <eigen3/Eigen/Dense>
#define M_PI 3.14159265358979323846  /* pi */

const int N  = 1000;
double x0    = 0.0;
double x1    = 100.0;
double k     = 1.0;

using namespace Eigen;
using namespace std;

int main() {

    // The eigen approach
    ArrayXd n                = ArrayXd::LinSpaced(N+1,0,N);
    double multiplier        = M_PI/N;
    Array<double, 1, N+1> nT = n.transpose();
    ArrayXd x                = - cos(multiplier*n);
    ArrayXd xsub             = x.middleRows(1, N-1);
    ArrayXd ysub             = (x1-x0)/2*xsub + (x1+x0)/2;

    ArrayXXd T               = cos((acos(x).matrix()*nT.matrix()).array());
    ArrayXXd Tsub            = cos((acos(xsub).matrix()*nT.matrix()).array());
    ArrayXd sqinx            = 1/sqrt(1-xsub*xsub);

    MatrixXd inv1x2          = (sqinx).matrix().asDiagonal();

    // Can't use the following to test elements of inv1x2
    // std::cout << inv1x2(0,0) << "\n";

    MatrixXd Usub            = inv1x2 * sin(((acos(xsub).matrix())*nT.matrix()).array()).matrix();
    MatrixXd dTsub           = Usub*(n.matrix().asDiagonal());
    MatrixXd d2Tsub          = ((sqinx*sqinx).matrix().asDiagonal())*((xsub.matrix().asDiagonal()) * (dTsub.matrix()) - (Tsub.matrix()) * ((n*n).matrix().asDiagonal()));

    MatrixXd d2T(N+1, N+1);
    RowVectorXd a            = (pow((-1),nT))*(nT*nT+1)*(nT*nT)/3;
    RowVectorXd b            = (nT*nT+1)*(nT*nT)/3;
    d2T.middleRows(1,N-1)    = d2Tsub; 
    d2T.row(0)               = a;
    d2T.row(N)               = b;

    MatrixXd D2              = d2T.matrix() * ((T.matrix()).inverse());
    MatrixXd E2              = D2.middleRows(1,N-1).middleCols(1,N-1);
    MatrixXd Y               = ysub.matrix().asDiagonal();
    MatrixXd H               = - (4 / ((x1-x0)*(x1-x0))) * E2 + k*Y;

    Eigen::EigenSolver<Eigen::MatrixXd> HE(H);
    VectorXcd D              = HE.eigenvalues();
    MatrixXcd V              = HE.eigenvectors();
    std::cout << HE.info() << std::endl;

    // Open ofstream
    ofstream Dfile;
    Dfile.open("D-output.txt");

    ofstream Vfile;
    Vfile.open("V-output.txt");

    ofstream V544file;
    V544file.open("V544-output.txt");

    Dfile.precision(15);
    Dfile << D.real() << "\n";

    Vfile.precision(15);
    Vfile << V.real() << "\n";

    V544file.precision(15);

	for(int i = 1; i<N-1; i++)
    {
		V544file << ysub[i-1];
        V544file << " "        << V.col(544).row(i-1).real() << "\n";
	}
    Dfile.close();
    Vfile.close();
	V544file.close();
	system("gnuplot -p plot.gp");
	system("rsvg-convert -w 2000 -o V544-plot.png V544-plot.svg");

}
