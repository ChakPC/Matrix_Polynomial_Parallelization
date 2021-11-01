// Serial Implementation of Schur Decomposition
#include "../complexMatrix.h"
using namespace std;

// Function to compute Givens matrix for given i, j
// Theta computed as such that after rotation, the value at (i, j) becomes 0
complexMatrix computeGivens(int i, int j, int sizeOfMatrix, complexNumber Aij, complexNumber Ajj) {
    // Compute theta
    complexNumber theta = atan((complexNumber(-1.0, 0.0) * Aij) / Ajj);
    complexNumber cosTheta = cos(theta);
    complexNumber sinTheta = sin(theta);

    // Compute Givens matrix
    complexMatrix Givens(sizeOfMatrix, vector<complexNumber>(sizeOfMatrix, complexNumber(0.0, 0.0)));
    for (int x = 0; x < sizeOfMatrix; x++) {
        Givens[x][x] = complexNumber(1.0, 0.0);
    }
    Givens[i][i] = cosTheta;
    Givens[j][j] = cosTheta;
    Givens[i][j] = sinTheta;
    Givens[j][i] = -1.0 * sinTheta;
    return Givens;
}

vector<complexMatrix> computeQR(complexMatrix& inputMatrix) {
    int sizeOfMatrix = inputMatrix.size();
    complexMatrix R(sizeOfMatrix, vector<complexNumber>(sizeOfMatrix, complexNumber(0.0, 0.0)));
    complexMatrix QTranspose = identityMatrix(sizeOfMatrix);

    // Perform rotations to form upper triangular matrix
    for (int i = 1; i < sizeOfMatrix; i++) {
        for (int j = 0; j < i; j++) {
            auto G = computeGivens(i, j, sizeOfMatrix, inputMatrix[i][j], inputMatrix[j][j]);
            // printf("\nG: %d, %d \n", i, j);
            // printMatrix(G);
            QTranspose = G * QTranspose;
        }
    }
    // printMatrix(QTranspose);

    // Compute Q and R
    complexMatrix Q = transposeMatrix(QTranspose);
    R = QTranspose * inputMatrix;
    return {Q, R};
}

// Function to perform iterations of Schur Decomposition
complexMatrix schurDecompositionSerial(complexMatrix& inputMatrix, int numberOfIterations) {
    for (int i = 0; i < numberOfIterations; i++) {
        auto QR_Vector = computeQR(inputMatrix);
        auto Q = QR_Vector[0];
        auto R = QR_Vector[1];
        // printf("\nQ\n");
        // printMatrix(Q);
        // printf("\nR\n");
        // printMatrix(R);
        inputMatrix = R * Q;
    }

    return inputMatrix;
}
// A = Qh * U * Q
// Q: unitary matrix
// Q * Qh = I
// Qh: conjugate transpose of Q
// QR method
// A0 = A
// A0 => Q0 * R0
// R0 * Q0 => A1
// A1 => Q1 * R1
// R1 * Q1 => A2
//....A = U
// Gram Schmidt
// ....
// Givens Rotations
/*
[[1 2 3]
 [4 5 6]
 [7 8 9]]

theta = atan(-Aij / Ajj)

G(i, j, theta)
=[[]
  []
  []]

*/

// Driver function for testing
int main() {
    // Take Input
    int n;
    cin >> n;
    complexMatrix A(n, vector<complexNumber>(n, {0.0, 0.0}));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    complexMatrix t = schurDecompositionSerial(A, 10);
    printMatrix(t);
    // complexMatrix t = computeGivens(2, 1, 3, complexNumber(8, 0), complexNumber(5, 0));
    // printMatrix(t);
    return 0;
}
