#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
using namespace std;

/**
 * Find the absolute value of a specified vector
 * 
 * @param vec      Vector whose absolute value is to be determined.
 * @return         abs(vec)
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
 * @param vvec     Vector of vectors that are to be added.
 * @return         Vector that is the sum of the vectors contained within vvec.
 */
vector<double> vecAdd(vector<vector<double>> vvec) {
    int N = vvec.size();
    int M = (vvec[0]).size();
    vector<double> addedVec(M);
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++) {
            addedVec[j] += vvec[i][j];
        }
    }
    return addedVec;
}

vector<double> vecAdd2(vector<double> vec1, vector<double> vec2) {
    int N = vec1.size();
    std::vector<double> vec3; 
    std::transform(vec1.begin(), vec1.end(), vec2.begin(), std::back_inserter(vec3), std::plus<float>());
    return vec3;
}

vector<double> vecAdd3(vector<double> vec1, vector<double> vec2, vector<double> vec3) {
    int N = vec1.size();
    std::vector<double> addedVec12 (N, 0.0);
    std::vector<double> addedVec123 (N, 0.0); 
    std::transform(vec1.begin(), vec1.end(), vec2.begin(), std::back_inserter(addedVec12), std::plus<float>());
    std::transform(addedVec12.begin(), addedVec12.end(), vec3.begin(), std::back_inserter(addedVec123), std::plus<float>());
    return addedVec123;
}

vector<double> vecAddTrans(vector<vector<double>> vvec) {
    int N = vvec.size();
    int M = vvec[0].size();
    vector<double> addedVec(M);
    for (int i = 0; i < N; i++) {
        std::transform(vvec[i].begin(), vvec[i].end(), addedVec.begin(), std::back_inserter(addedVec), std::plus<float>());
    }
    return addedVec;
}

/**
 * Find the maximum value of a specified vector
 *
 * @param vec      Vector whose maximum value is being found.
 * @return         Maximum element of vec.
 */
double vecMax(vector<double> vec) {
    double max = *std::max_element(vec.begin(), vec.end());
    return max;
}

/**
 * Find the minimum value of a specified vector
 *
 * @param vec      Vector whose minimum value is being found.
 * @return         Maximum element of vec.
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
    // if (N1 != N2) {
    //     throw std::invalid_argument("Argument vectors must be of the same size!");
    // }

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
// vector<double> vecScalMult(vector<double> vec, double scalar) {
//     using namespace std::placeholders;
//     vector<double> vec2 = for_each(vec.begin(), vec.end(), [scalar](double &c){ c *= scalar; });
//     return vec2;
// }