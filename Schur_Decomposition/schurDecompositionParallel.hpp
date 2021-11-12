#include "../complexMatrixParallel.hpp"
/*
 * Function to compute Givens matrix for given i, j
 * Theta computed as such that after rotation, the value at (i, j) becomes 0
 * @param i row (0 based indexing)
 * @param j column (0 based indexing)
 * @param sizeOfMatrix size of matrix on which givens is to be calculated
 * @param Aij value in i'th row and j'th column of Matrix
 * @param Ajj value in j'th row and j'th column of Matrix
 * @return Given's matrix
 */
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

/*
 * multiplication with Given's Matrix
 * @param G Given's Matrix
 * @param operand Matrix to multiply with Given's Matrix
 * @param r row number (0 based)
 * @param c column number (0 based)
 */
inline complexMatrix multG(complexMatrix &G, complexMatrix &operand, int r, int c) {
    complexMatrix res = operand;
    complexNumber sumR = complexNumber(0, 0);
    complexNumber sumC = complexNumber(0, 0);
    int cols = operand[0].size();
#pragma omp parallel for
    for (int i = 0; i < cols; i++) {
        res[r][i] = G[r][c] * operand[c][i] + G[r][r] * operand[r][i];
        res[c][i] = G[c][c] * operand[c][i] + G[c][r] * operand[r][i];
    }
    return res;
}

/*
 * compute QR decomposition of matrix
 * @param inputMatrix matrix for which QR decomposition is to be calculated
 * @return vector of complex Number, [0] => Q Matrix, [1] => R Matrix
 */
inline vector<complexMatrix> computeQRParallel(complexMatrix &inputMatrix) {
    int sizeOfMatrix = inputMatrix.size();
    complexMatrix R = inputMatrix;
    complexMatrix QTranspose = identityMatrixParallel(sizeOfMatrix);
    for (int diff = 1; diff < sizeOfMatrix; diff++) {
        int limit = sizeOfMatrix - diff - 1;
        int flag = 0;
#pragma omp parallel for
        for (int col = 0; col < limit; col += 2) {
            int row = diff + col;
            if (inputMatrix[row][col] != complexNumber(0, 0)) {
                complexNumber factor = R[row][col] / R[col][col];
                if (factor == complexNumber(0, 1) || factor == complexNumber(0, -1)) {
                    flag = 1;
                }
                if (!flag) {
                    auto G = computeGivensParallel(row, col, sizeOfMatrix, R[row][col], R[col][col]);
#pragma omp critical
                    {
                        QTranspose = multG(G, QTranspose, row, col);
                        R = multG(G, R, row, col);
                    }
                }
            }
        }
        if (flag) {
            return {};
        }
        flag = 0;
#pragma omp parallel for
        for (int col = 1; col < limit; col += 2) {
            int row = diff + col;
            if (inputMatrix[row][col] != complexNumber(0, 0)) {
                complexNumber factor = R[row][col] / R[col][col];
                if (factor == complexNumber(0, 1) || factor == complexNumber(0, -1)) {
                    flag = 1;
                }
                if (!flag) {
                    auto G = computeGivensParallel(row, col, sizeOfMatrix, R[row][col], R[col][col]);
#pragma omp critical
                    {
                        QTranspose = multG(G, QTranspose, row, col);
                        R = multG(G, R, row, col);
                    }
                }
            }
        }
        if (flag) {
            return {};
        }
    }

    // Compute Q and R
    complexMatrix Q = transposeMatrixParallel(QTranspose);
    return {Q, R};
}

/*
 * Function to perform iterations of Schur Decomposition
 * @param inputMatrix matrix for calculating Schur Decomposition
 * @param numberOfIterations total iterations to perform schur decomposition
 * @return vector of complex Matrix, [0] => Q (unitary matrix), [1] => upper triangular matrix
 */
inline vector<complexMatrix> schurDecompositionParallel(complexMatrix inputMatrix, int numberOfIterations) {
    complexMatrix Q, R;
    complexMatrix QFinal = identityMatrixParallel(inputMatrix.size());
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
    processZeroParallel(QFinal);
    processZeroParallel(inputMatrix);
    return {transposeMatrixParallel(QFinal), inputMatrix};
}
