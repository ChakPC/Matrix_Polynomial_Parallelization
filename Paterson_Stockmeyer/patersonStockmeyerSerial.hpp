// Header file for serial Implementation of Paterson - Stockmeyer algorithm

// Function to compute powers of input matrix till exponent limit
// Brute Force Implmentation: O(limit * n * n)
inline vector<complexMatrix> computePowersSerial(complexMatrix &inputMatrix, int limit) {
    int sz = inputMatrix.size();
    vector<complexMatrix> res(limit + 1, complexMatrix(sz, vector<complexNumber>(sz, complexNumber(0.0, 0.0))));
    complexMatrix temp = inputMatrix;
    complexMatrix identity = identityMatrixSerial(sz);
    res[0] = identity;
    res[1] = inputMatrix;
    for (int i = 2; i <= limit; i++) {
        temp = multiplySerial(temp, inputMatrix);
        res[i] = temp;
    }
    return res;
}

// Serial Implementation of Paterson Stockmeyer
inline complexMatrix patersonStockmeyerSerial(complexMatrix &inputMatrix, vector<complexNumber> &coefficients, int polynomialVariable, int polynomialDegree) {
    int inputMatrixSize = inputMatrix.size();
    int degree = coefficients.size() - 1;

    // Compute powers of A till p
    vector<complexMatrix> powersOfInputMatrix = computePowersSerial(inputMatrix, polynomialVariable);

    // Declare resultant matrix
    complexMatrix resultantMatrix(inputMatrixSize, vector<complexNumber>(inputMatrixSize, complexNumber(0.0, 0.0)));

    // Horner's loop
    for (int q = polynomialDegree - 1; q >= 0; q--) {
        complexMatrix temp(inputMatrixSize, vector<complexNumber>(inputMatrixSize, complexNumber(0.0, 0.0)));
        for (int j = 0; j < polynomialVariable; j++) {
            complexNumber coeff = {0.0, 0.0};
            int index = polynomialVariable * q + j;
            if (j <= degree) coeff = coefficients[index];
            complexMatrix adder = multiplySerial(coeff, powersOfInputMatrix[j]);
            temp = addSerial(temp, adder);
        }
        complexMatrix adder = multiplySerial(resultantMatrix, powersOfInputMatrix[polynomialVariable]);
        resultantMatrix = addSerial(adder, temp);
    }

    // Return answer
    processZeroSerial(resultantMatrix);
    return resultantMatrix;
}