// Header file for serial Implementation of Schur Decomposition

// Function to compute Givens matrix for given i, j
// Theta computed as such that after rotation, the value at (i, j) becomes 0
inline complexMatrix computeGivensParallel(int i, int j, int sizeOfMatrix, complexNumber Aij, complexNumber Ajj) {
    // Compute theta
    complexNumber theta = atan((complexNumber(-1.0, 0.0) * Aij) / (Ajj));
    complexNumber cosTheta = cos(theta);
    complexNumber sinTheta = sin(theta);

    // Compute Givens matrix
    complexMatrix Givens(sizeOfMatrix, vector<complexNumber>(sizeOfMatrix, complexNumber(0.0, 0.0)));

#pragma omp parallel for
    for (int x = 0; x < sizeOfMatrix; x++) {
        Givens[x][x] = complexNumber(1.0, 0.0);
    }

    Givens[i][i] = cosTheta;
    Givens[j][j] = cosTheta;
    Givens[i][j] = sinTheta;
    Givens[j][i] = -1.0 * sinTheta;
    return Givens;
}

inline vector<complexMatrix> computeQRParallel(complexMatrix &inputMatrix) {
    int sizeOfMatrix = inputMatrix.size();
    complexMatrix R = inputMatrix;
    complexMatrix QTranspose = identityMatrix(sizeOfMatrix);
    for (int i = 1; i < sizeOfMatrix; i++) {
        for (int j = 0; j < i; j++) {
            if (inputMatrix[i][j] == complexNumber(0, 0)) continue;
            complexNumber factor = R[i][j] / R[j][j];
            if (factor == complexNumber(0, 1) || factor == complexNumber(0, -1)) {
                return {};
            }
            auto G = computeGivensParallel(i, j, sizeOfMatrix, R[i][j], R[j][j]);
            QTranspose = G * QTranspose;
            R = G * R;
        }
    }

    // Compute Q and R
    complexMatrix Q = transposeMatrix(QTranspose);
    return {Q, R};
}

// Function to perform iterations of Schur Decomposition
inline vector<complexMatrix> schurDecompositionParallel(complexMatrix inputMatrix, int numberOfIterations) {
    complexMatrix Q, R;
    complexMatrix QFinal = identityMatrix(inputMatrix.size());
    for (int i = 0; i < numberOfIterations; i++) {
        vector<complexMatrix> QR_Vector = computeQRParallel(inputMatrix);
        if (QR_Vector.size() == 0) {
            return {};
        }
        Q = QR_Vector[0];
        R = QR_Vector[1];
        inputMatrix = R * Q;
        QFinal = QFinal * Q;
    }
    processZeroSerial(QFinal);
    processZeroSerial(inputMatrix);
    return {transposeMatrix(QFinal), inputMatrix};
}
