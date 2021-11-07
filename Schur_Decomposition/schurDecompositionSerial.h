// Header file for serial Implementation of Schur Decomposition

// Function to compute Givens matrix for given i, j
// Theta computed as such that after rotation, the value at (i, j) becomes 0
inline complexMatrix computeGivens(int i, int j, int sizeOfMatrix, complexNumber Aij, complexNumber Ajj)
{
    // Compute theta
    complexNumber theta = atan((complexNumber(-1.0, 0.0) * Aij) / Ajj);
    complexNumber cosTheta = cos(theta);
    complexNumber sinTheta = sin(theta);

    // Compute Givens matrix
    complexMatrix Givens(sizeOfMatrix, vector<complexNumber>(sizeOfMatrix, complexNumber(0.0, 0.0)));
    for (int x = 0; x < sizeOfMatrix; x++)
    {
        Givens[x][x] = complexNumber(1.0, 0.0);
    }
    Givens[i][i] = cosTheta;
    Givens[j][j] = cosTheta;
    Givens[i][j] = sinTheta;
    Givens[j][i] = -1.0 * sinTheta;
    return Givens;
}

inline vector<complexMatrix> computeQR(complexMatrix &inputMatrix)
{
    int sizeOfMatrix = inputMatrix.size();
    complexMatrix R=inputMatrix;
    complexMatrix QTranspose = identityMatrix(sizeOfMatrix);
    for (int i = 1; i < sizeOfMatrix; i++)
    {
        for (int j = 0; j < i; j++)
        {
            auto G = computeGivens(i, j, sizeOfMatrix, R[i][j], R[j][j]);
            QTranspose = G * QTranspose;
            R = G*R;
        }
    }

    // Compute Q and R
    complexMatrix Q = transposeMatrix(QTranspose);
    return {Q, R};
}

// Function to perform iterations of Schur Decomposition
inline vector<complexMatrix> schurDecompositionSerial(complexMatrix &inputMatrix, int numberOfIterations)
{
    complexMatrix Q, R;
    complexMatrix QFinal = identityMatrix(inputMatrix.size());
    for (int i = 0; i < numberOfIterations; i++)
    {
        vector<complexMatrix> QR_Vector = computeQR(inputMatrix);
        Q = QR_Vector[0];
        R = QR_Vector[1];
        inputMatrix = R * Q;
        QFinal = QFinal * Q;
    }
    return {transposeMatrix(QFinal), inputMatrix};
}

