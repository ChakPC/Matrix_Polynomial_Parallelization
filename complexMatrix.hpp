// Header file for complex matrix operations

#include <bits/stdc++.h>
#include <omp.h>

using namespace std;

#define complexNumber complex<double>
#define complexMatrix vector<vector<complexNumber>>
#define clusterThreshold 0.1
#define zeroLimit 1e-5

bool areSame(complexNumber x, complexNumber y) {
    return (fabs(real(x) - real(y)) < zeroLimit) && (fabs(imag(x) - imag(y)) < zeroLimit);
}

// + operator overloaded to add two complex matrices
inline complexMatrix operator+(complexMatrix &A, complexMatrix &B) {
    int row = A.size(), col = A[0].size();
    complexMatrix res(row, vector<complexNumber>(col, complexNumber(0.0, 0.0)));

#pragma omp parallel for collapse(2)
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

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

// * operator overloaded to multiply a complex number with a complex matrix
inline complexMatrix operator*(complexNumber coeff, complexMatrix &matrix) {
    int row = matrix.size(), col = matrix[0].size();
    complexMatrix res(row, vector<complexNumber>(col, {0.0, 0.0}));

#pragma omp parallel for collapse(2)
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            res[i][j] = coeff * matrix[i][j];
        }
    }
    return res;
}

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

// * operator overloaded to multiply two complex matrices
inline complexMatrix operator*(complexMatrix &A, complexMatrix &B) {
    int rowA = A.size(), colA = A[0].size(), rowB = B.size(), colB = B[0].size();
    complexMatrix res(rowA, vector<complexNumber>(colB, {0.0, 0.0}));

#pragma omp parallel for collapse(2)
    for (int i = 0; i < rowA; i++) {
        for (int j = 0; j < colB; j++) {
            for (int k = 0; k < colA; k++) {
                res[i][j] = res[i][j] + A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

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

inline void processZeroParallel(complexMatrix &inputMatrix) {
    int row = inputMatrix.size(), col = inputMatrix[0].size();

#pragma omp parallel for collapse(2)
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

bool areMatricesEqual(complexMatrix &A, complexMatrix &B) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) return false;
    volatile bool flag = false;
#pragma omp parallel for shared(flag) collapse(2)
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

// Function to print complex matrix
inline void printMatrix(complexMatrix &res) {
    int n = res.size();
    if (n == 0) return;
    int m = res[0].size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << "( " << real(res[i][j]) << " ) + i( " << imag(res[i][j]) << " )  ";
        }
        cout << "\n";
    }
}

// Function to print complex number
inline void printComplex(complexNumber res) {
    cout << "(" << real(res) << " + i" << imag(res) << ") ";
}

// Function to return transpose of a complex matrix
// Returns the conjugate transpose
inline complexMatrix transposeMatrix(complexMatrix &inputMatrix) {
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

inline complexMatrix identityMatrix(int n) {
    complexMatrix res(n, vector<complexNumber>(n, complexNumber(0.0, 0.0)));
    for (int i = 0; i < n; i++) {
        res[i][i] = complexNumber(1.0, 0.0);
    }
    return res;
}

// Function to compute kronecker product of A and B
inline complexMatrix computeKroneckerProduct(complexMatrix &A, complexMatrix &B) {
    int rA = A.size(), cA = A[0].size();
    int rB = B.size(), cB = B[0].size();
    complexMatrix res(rA * rB, vector<complexNumber>(cA * cB, complexNumber(0, 0)));
    for (int i = 0; i < rA; i++) {
        for (int j = 0; j < cA; j++) {
            complexMatrix Aij_B = A[i][j] * B;
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

// Function to stack columns of a matrix atop each other from left to right
inline complexMatrix columnStack(complexMatrix &inputMatrix) {
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

// Function to convert column stacked matrix back to standard form
inline complexMatrix reverseStacking(complexMatrix &inputMatrix, int rows) {
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

// Comparator to compare complex number
inline bool cmp(complexNumber &A, complexNumber &B) {
    double a = real(A) * real(A) + imag(A) * imag(A);
    double b = real(B) * real(B) + imag(B) * imag(B);
    return a < b;
}

// Check for clustered Eigen Values
inline bool checkClusteredEigenValues(complexMatrix &T) {
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

// Function to update a row for gauss elimination
inline void updateRow(complexMatrix &A, complexMatrix &B, int i, int j) {
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

// Function to perform Gauss Elimination
inline complexMatrix gaussElimination(complexMatrix &A, complexMatrix &B) {
    // Ax = B
    int size = A.size();
    for (int i = 0; i < size - 1; i++) {
        for (int j = i + 1; j < size; j++) {
            updateRow(A, B, i, j);
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