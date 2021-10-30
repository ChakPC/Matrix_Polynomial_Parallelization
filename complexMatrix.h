// Header file for complex matrix operations

#include <bits/stdc++.h>
using namespace std;

#define complexNumber complex<double>
#define complexMatrix vector<vector<complexNumber>>

// + operaotr overloaded to add two complex matrices
inline complexMatrix operator+(complexMatrix &A, complexMatrix &B){
    int size = A.size();
    complexMatrix res(size, vector<complexNumber>(size, complexNumber(0.0, 0.0)));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

// * operator overloaded to multiply a complex number with a complex matrix
inline complexMatrix operator*(complexNumber coeff, complexMatrix &matrix){
    int size = matrix.size();
    complexMatrix res(size, vector<complexNumber>(size, {0.0, 0.0}));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            res[i][j] = coeff * matrix[i][j];
        }
    }
    return res;
}

// * operator overloaded to multiply two complex matrices
inline complexMatrix operator*(complexMatrix &A, complexMatrix &B){
    int size = A.size();
    complexMatrix res(size, vector<complexNumber>(size, {0.0, 0.0}));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                res[i][j] = res[i][j] + A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

// Function to print complex matrix
inline void printMatrix(complexMatrix &res){
    int n = res.size();
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << "( " << real(res[i][j]) <<" ) + i( "<< imag(res[i][j]) <<" )  ";
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
inline complexMatrix transposeMatrix(complexMatrix& inputMatrix) {
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
    for(int i=0; i<n; i++){
        res[i][i] = complexNumber(1.0, 0.0);
    }
    return res;
}
