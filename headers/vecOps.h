#include <vector>
#include <cmath>
using namespace std;

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
 * Find the minimum value of a specified vector
 *
 * @param vec           Vector whose minimum value is being found.
 * @return              Maximum element of vec.
 */
double vecMin(vector<double> vec) {
    double min = *std::min_element(vec.begin(), vec.end());
    return min;
}

/**
 * Element-wise multiplication between two vectors
 * 
 * @param vec1     First vector to be multiplied.
 * @param vec2     Second vector to be multiplied.
 * @return         vec1 * vec2 (element-wise)
 */
vector<double> vecMult(vector<double> vec1, vector<double> vec2) {
    // Initialize relevant variables
    int N1 = vec1.size();
    int N2 = vec2.size();
    vector<double> multvec (N1, 0.0);

    // Return an error if N1 != N2
    if (N1 != N2) {
        throw std::invalid_argument("Argument vectors must be of the same size!");
    }

    // Loop over each element, multiply them and write the result to multvec
    for (int i = 0; i < N1; i++) {
        multvec[i] = vec1[i]*vec2[i];
    }
    return multvec;
}

/**
 * Multiply a vector by a scalar
 * 
 * @param vec      Vector to be multiplied.
 * @param scalar   A floating point number the vector is to be multiplied by.
 * @return         vec * scalar
 */
vector<double> vecScalMult(vector<double> vec, double scalar) {
    // Initialize relevant variables
    int N = vec.size();
    vector<double> multvec (N, 0.0);

    // Loop over each vector element and multiply it by the scalar
    for (int i = 0; i < N; i++) {
        multvec[i] = vec[i]*scalar;
    }
    return multvec;
}