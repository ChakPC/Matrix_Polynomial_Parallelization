#include <bits/stdc++.h>

#define complexNumber complex<double>
using namespace std;

vector<vector<complexNumber>> operator*(vector<vector<complexNumber>>& A, vector<vector<complexNumber>>& B) {
    int size = A.size();
    vector<vector<complexNumber>> res(size, vector<complexNumber>(size, {0.0, 0.0}));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            for (int k = 0; k < size; k++) {
                res[i][j] = res[i][j] + A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

vector<vector<complexNumber>> transposeMatrix(vector<vector<complexNumber>>& inputMatrix) {
    int n = inputMatrix.size();
    vector<vector<complexNumber>> newMatrix(n, vector<complexNumber>(n, complexNumber(0, 0)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            newMatrix[i][j] = inputMatrix[j][i];
        }
    }
    return newMatrix;
}

void printComplex(complexNumber res) {
    cout << "(" << real(res) << " + i" << imag(res) << ") ";
}

void printMatrix(vector<vector<complexNumber>> res) {
    int n = res.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "(" << real(res[i][j]) << " + i" << imag(res[i][j]) << ") ";
        }
        cout << "\n";
    }
}
void identityMatrix(vector<vector<complexNumber>>& matrix) {
    int n = matrix.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i != j)
                matrix[i][j] = complexNumber(0.0, 0);
            else
                matrix[i][j] = complexNumber(1.0, 0);
        }
    }
}
vector<vector<complexNumber>> computeGivens(int i, int j, int sizeOfMatrix, complexNumber Aij, complexNumber Ajj) {
    complexNumber theta = atan(complexNumber(-1.0, 0.0) * Aij / Ajj);
    auto cosTheta = cos(theta);
    auto sinTheta = sin(theta);
    vector<vector<complexNumber>> Givens(sizeOfMatrix, vector<complexNumber>(sizeOfMatrix, complexNumber(0.0, 0.0)));
    for (int x = 0; x < sizeOfMatrix; x++) {
        Givens[x][x] = complexNumber(1.0, 0.0);
    }
    Givens[i][i] = cosTheta;
    Givens[j][j] = cosTheta;
    Givens[i][j] = sinTheta;
    Givens[j][i] = -sinTheta;
    return Givens;
}

vector<vector<vector<complexNumber>>> computeQR(vector<vector<complexNumber>>& inputMatrix) {
    int sizeOfMatrix = inputMatrix.size();
    vector<vector<complexNumber>> Q(sizeOfMatrix, vector<complexNumber>(sizeOfMatrix, complexNumber(0.0, 0.0)));
    vector<vector<complexNumber>> R(sizeOfMatrix, vector<complexNumber>(sizeOfMatrix, complexNumber(0.0, 0.0)));
    identityMatrix(Q);

    for (int i = 1; i < sizeOfMatrix; i++) {
        for (int j = 0; j < i; j++) {
            auto G = computeGivens(i, j, sizeOfMatrix, inputMatrix[i][j], inputMatrix[j][j]);
            Q = Q * G;
            printMatrix(G * inputMatrix);
        }
    }
    auto transpose = transposeMatrix(Q);
    R = transpose * inputMatrix;
    return {Q, R};
}

vector<vector<complexNumber>> schurDecomposition(vector<vector<complexNumber>>& inputMatrix, int numberOfIterations) {
    for (int i = 0; i < numberOfIterations; i++) {
        auto QR_Vector = computeQR(inputMatrix);
        auto Q = QR_Vector[0];
        auto R = QR_Vector[1];
        // cout << endl;
        // printMatrix(Q);
        // printMatrix(R);
        // cout << endl;
        inputMatrix = R * Q;
    }
    return inputMatrix;
}

int main() {
    vector<vector<complexNumber>> inputMatrix = {
        {complexNumber(12, 0), complexNumber(-51, 0), complexNumber(4, 0)},
        {complexNumber(6, 0), complexNumber(167, 0), complexNumber(-68, 0)},
        {complexNumber(-4, 0), complexNumber(24, 0), complexNumber(-41, 0)}};

    printMatrix(inputMatrix);
    printMatrix(schurDecomposition(inputMatrix, 1));
    return 0;
}
