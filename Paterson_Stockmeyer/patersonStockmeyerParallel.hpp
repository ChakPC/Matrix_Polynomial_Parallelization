// Header file for serial Implementation of Paterson - Stockmeyer algorithm

// Function to compute powers of input matrix till exponent limit
// Brute Force Implmentation: O(limit * n * n)
inline vector<complexMatrix> computePowersParallel(complexMatrix &inputMatrix, int limit) {
    int sz = inputMatrix.size();
    vector<complexMatrix> res(limit + 1, complexMatrix(sz, vector<complexNumber>(sz, complexNumber(0.0, 0.0))));
    complexMatrix temp = inputMatrix;
    complexMatrix identity(sz, vector<complexNumber>(sz, complexNumber(0.0, 0.0)));
    for (int i = 0; i < sz; i++) {
        identity[i][i] = complexNumber(1.0, 0.0);
    }
    res[0] = identity;
    res[1] = inputMatrix;
    for (int i = 2; i <= limit; i++) {
        temp = temp * inputMatrix;
        res[i] = temp;
    }
    return res;
}

// Parallel Implementation of Paterson Stockmeyer
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

    // Return answer
    processZeroParallel(resultantMatrix);
    return resultantMatrix;
}