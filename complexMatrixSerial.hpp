/*
 * add complex matrices serially
 * @param A Matrix A
 * @param B Matrix B
 * @return A+B
 */
complexMatrix addSerial(complexMatrix &A, complexMatrix &B) {
    int row = A.size(), col = A[0].size();
    complexMatrix res(row, vector<complexNumber>(col, complexNumber(0.0, 0.0)));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

/*
 * multiply complex matrices serially
 * @param A coefficient A
 * @param B Matrix B
 * @return A*B
 */
complexMatrix multiplySerial(complexNumber coeff, complexMatrix &matrix) {
    int row = matrix.size(), col = matrix[0].size();
    complexMatrix res(row, vector<complexNumber>(col, {0.0, 0.0}));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            res[i][j] = coeff * matrix[i][j];
        }
    }
    return res;
}

/*
 * multiply complex matrices serially
 * @param A Matrix A
 * @param B Matrix B
 * @return A*B
 */
complexMatrix multiplySerial(complexMatrix &A, complexMatrix &B) {
    int rowA = A.size(), colA = A[0].size(), rowB = B.size(), colB = B[0].size();
    complexMatrix res(rowA, vector<complexNumber>(colB, {0.0, 0.0}));
    for (int i = 0; i < rowA; i++) {
        for (int j = 0; j < colB; j++) {
            for (int k = 0; k < colA; k++) {
                res[i][j] = res[i][j] + A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

/*
 * Transforms matrix values to zero if they are smaller than zeroLimit
 * @param inputMatrix The input matrix to transform
 */
inline void processZeroSerial(complexMatrix &inputMatrix) {
    int row = inputMatrix.size(), col = inputMatrix[0].size();
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            if (abs(real(inputMatrix[i][j])) <= zeroLimit) {
                inputMatrix[i][j] = complexNumber(0, imag(inputMatrix[i][j]));
            }
            if (abs(imag(inputMatrix[i][j])) <= zeroLimit) {
                inputMatrix[i][j] = complexNumber(real(inputMatrix[i][j]), 0);
            }
        }
    }
}

/*
 * check if two matrices are equal by checking if the difference between two
 * values is smaller than zeroLimit
 * @param A Matrix A
 * @param B Matrix B
 * @return true if A == B else false
 */
bool areMatricesEqualSerial(complexMatrix &A, complexMatrix &B) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) return false;
    bool flag = false;
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            if (flag) continue;
            if (!areSame(A[i][j], B[i][j])) {
                flag = true;
            }
        }
    }
    if (flag) return false;
    return true;
}

/*
 * calculate conjugate transpose of a matrix
 * @param inputMatrix matrix for which conjugate transpose is to be calculated
 * @return conjugate transpose of inputMatrix
 */
inline complexMatrix transposeMatrixSerial(complexMatrix &inputMatrix) {
    int n = inputMatrix.size();
    complexMatrix newMatrix(n, vector<complexNumber>(n, complexNumber(0, 0)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double Re = real(inputMatrix[j][i]);
            double Im = imag(inputMatrix[j][i]);
            newMatrix[i][j] = complexNumber(Re, -1.0 * Im);
        }
    }
    return newMatrix;
}

/*
 * create a new identity matrix
 * @param n size of identity matrix
 * @return new identity matrix of size n
 */
inline complexMatrix identityMatrixSerial(int n) {
    complexMatrix res(n, vector<complexNumber>(n, complexNumber(0.0, 0.0)));
    for (int i = 0; i < n; i++) {
        res[i][i] = complexNumber(1.0, 0.0);
    }
    return res;
}

/*
 * calculate Kronecker Product of matrices
 * @param A Matrix A
 * @param B Matrix B
 * @return Kronecker Product of A and B
 */
inline complexMatrix computeKroneckerProductSerial(complexMatrix &A, complexMatrix &B) {
    int rA = A.size(), cA = A[0].size();
    int rB = B.size(), cB = B[0].size();
    complexMatrix res(rA * rB, vector<complexNumber>(cA * cB, complexNumber(0, 0)));
    for (int i = 0; i < rA; i++) {
        for (int j = 0; j < cA; j++) {
            complexMatrix Aij_B = multiplySerial(A[i][j], B);
            int rowStart = i * rB, rowEnd = (i + 1) * rB;
            int colStart = j * cB, colEnd = (j + 1) * cB;
            for (int x = rowStart; x < rowEnd; x++) {
                for (int y = colStart; y < colEnd; y++) {
                    res[x][y] = Aij_B[x - rowStart][y - colStart];
                }
            }
        }
    }
    processZeroSerial(res);
    return res;
}

/*
 * stack columns of a matrix atop each other from left to right
 * @param inputMatrix matrix for stacking
 * @return column Stacked Matrix of inputMatrix
 */
inline complexMatrix columnStackSerial(complexMatrix &inputMatrix) {
    int p = inputMatrix.size(), q = inputMatrix[0].size();
    complexMatrix res(p * q, vector<complexNumber>(1, complexNumber(0, 0)));
    int pos = p * q - 1;
    for (int i = 0; i < q; i++) {
        for (int j = p - 1; j >= 0; j--) {
            res[pos][0] = inputMatrix[j][i];
            pos--;
        }
    }
    return res;
}

/*
 * Function to convert column stacked matrix back to standard form
 * @param inputMatrix matrix for reverse stacking
 * @return reverse stacked matrix for inputMatrix
 */
inline complexMatrix reverseStackingSerial(complexMatrix &inputMatrix, int rows) {
    int prod = inputMatrix.size();
    int cols = prod / rows;
    complexMatrix res(rows, vector<complexNumber>(cols, complexNumber(0, 0)));
    int pos = 0;
    for (int i = cols - 1; i >= 0; i--) {
        for (int j = 0; j < rows; j++) {
            res[j][i] = inputMatrix[pos][0];
            pos++;
        }
    }
    return res;
}

/*
 * check if clustered eigen values are present in a matrix
 * @param T matrix
 * @return true if clustered eigen values are present else false
 */
inline bool checkClusteredEigenValuesSerial(complexMatrix &T) {
    int size = T.size();
    vector<complexNumber> eigenValues;
    for (int i = 0; i < size; i++) {
        if (i + 1 < size && abs(real(T[i + 1][i])) > zeroLimit) {
            double a = real(T[i][i]), b = real(T[i][i + 1]), c = real(T[i + 1][i]);
            eigenValues.push_back(complexNumber(a, sqrt(-b * c)));
            eigenValues.push_back(complexNumber(a, -sqrt(-b * c)));
            i++;
        } else {
            eigenValues.push_back(T[i][i]);
        }
    }
    sort(eigenValues.begin(), eigenValues.end(), cmp);
    for (int i = 0; i < eigenValues.size() - 1; i++) {
        double t1 = real(eigenValues[i + 1]) * real(eigenValues[i + 1]) + imag(eigenValues[i + 1]) * imag(eigenValues[i + 1]);
        double t2 = real(eigenValues[i]) * real(eigenValues[i]) + imag(eigenValues[i]) * imag(eigenValues[i]);
        if (t1 - t2 <= clusterThreshold) return 1;
    }
    return 0;
}

/*
 * update a row for gauss elimination
 * @param A Matrix A
 * @param B Matrix B
 * @param i row number (0 based)
 * @param j column number (0 based)
 */
inline void updateRowSerial(complexMatrix &A, complexMatrix &B, int i, int j) {
    int size = A.size();
    if (A[i][i] == complexNumber(0, 0)) {
        int candidateRow = -1;
        for (int k = i + 1; k < size; k++) {
            if (A[k][i] != complexNumber(0, 0)) {
                candidateRow = k;
                break;
            }
        }
        if (candidateRow == -1) {
            cout << "Sylvester Equation Not solvable :(";
            exit(0);
        }
        swap(A[i], A[candidateRow]);
        swap(B[i], B[candidateRow]);
    }
    complexNumber factor = A[j][i] / A[i][i];
    for (int col = 0; col < size; col++) {
        A[j][col] = A[j][col] - (factor * A[i][col]);
    }
    B[j][0] = B[j][0] - (factor * B[i][0]);
}

/*
 * perform Gauss Elimination
 * @param A Matrix A
 * @param B Matrix B
 * @return Gauss Elimination result
 */
inline complexMatrix gaussEliminationSerial(complexMatrix &A, complexMatrix &B) {
    // Ax = B
    int size = A.size();
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            updateRowSerial(A, B, i, j);
        }
    }
    processZeroSerial(A);
    complexMatrix res(size, vector<complexNumber>(1));
    for (int x = size - 1; x >= 0; x--) {
        complexNumber sum(0, 0);
        for (int y = x + 1; y < size; y++) {
            sum += A[x][y] * res[y][0];
        }
        res[x][0] = (B[x][0] - sum) / A[x][x];
    }
    processZeroSerial(res);
    return res;
}