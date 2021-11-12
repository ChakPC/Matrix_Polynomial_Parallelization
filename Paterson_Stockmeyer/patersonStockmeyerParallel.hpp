#include "../complexMatrixParallel.hpp"
/*
 * Function to compute powers of input matrix till exponent limit
 * @param inputMatrix matrix for calculating powers
 * @param limit exponent limit
 */
inline vector<complexMatrix> computePowersParallel(complexMatrix &inputMatrix, int limit) {
    int sz = inputMatrix.size();
    vector<complexMatrix> res(limit + 1, complexMatrix(sz, vector<complexNumber>(sz, complexNumber(0.0, 0.0))));
    complexMatrix temp = inputMatrix;
    complexMatrix identity = identityMatrixParallel(sz);
    map<int, int> vis;
    res[0] = identity;
    vis[0] = 1;
    res[1] = inputMatrix;
    vis[1] = 1;
    int logLimit = floor(log2(limit));
    int power = 2;
    for (int i = 1; i <= logLimit; i++) {
        temp = temp * temp;
        res[1 << i] = temp;
        vis[1 << i] = 1;
    }
#pragma omp parallel for
    for (int i = 2; i <= limit; i += 2) {
        if (!vis[i]) {
            vector<int> decomposition;
            int val = i;
            while (val > 0) {
                int greatestPower = floor(log2(val));
                decomposition.push_back(greatestPower);
                val -= 1 << greatestPower;
            }
            complexMatrix matrix = identity;
            for (int j = 0; j < decomposition.size(); j++) {
                matrix = matrix * res[1 << decomposition[j]];
            }
            res[i] = matrix;
            vis[i] = 1;
            if (i + 1 <= limit) {
                res[i + 1] = matrix * inputMatrix;
                vis[i + 1] = 1;
            }
        }
    }

    return res;
}

/*
 * Parallel Implementation of Paterson Stockmeyer
 * @param inputMatrix complex Matrix
 * @param coefficients vector of complex numbers
 * @param polynomialVariable exponent limit for computing powers
 * @param polynomialDegree degree of polynomial
 */
inline complexMatrix patersonStockmeyerParallel(complexMatrix &inputMatrix, vector<complexNumber> &coefficients, int polynomialVariable, int polynomialDegree) {
    int inputMatrixSize = inputMatrix.size();
    int degree = coefficients.size() - 1;

    // Compute powers of A till p
    vector<complexMatrix> powersOfInputMatrix = computePowersParallel(inputMatrix, polynomialVariable);

    // Declare resultant matrix
    complexMatrix resultantMatrix(inputMatrixSize, vector<complexNumber>(inputMatrixSize, complexNumber(0.0, 0.0)));

    // Horner's loop
    for (int q = polynomialDegree - 1; q >= 0; q--) {
        complexMatrix temp(inputMatrixSize, vector<complexNumber>(inputMatrixSize, complexNumber(0.0, 0.0)));

#pragma omp parallel for
        for (int j = 0; j < polynomialVariable; j++) {
            complexNumber coeff = {0.0, 0.0};
            int index = polynomialVariable * q + j;
            if (j <= degree) coeff = coefficients[index];
            complexMatrix adder = coeff * powersOfInputMatrix[j];

#pragma omp critical
            {
                temp = temp + adder;
            }
        }

        complexMatrix adder = resultantMatrix * powersOfInputMatrix[polynomialVariable];
        resultantMatrix = adder + temp;
    }

    processZeroParallel(resultantMatrix);
    return resultantMatrix;
}