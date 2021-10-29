// Serial Implementation of Paterson - Stockmeyer algorithm
#include <bits/stdc++.h>
#include <complex>
using namespace std;

struct complexNumber {
    double Re;
    double Im;
    complexNumber(double r = 0.0, double i = 0.0) {
        Re = r;
        Im = i;
    }
    complexNumber operator+(complexNumber operand) {
        complexNumber res;
        res.Re = Re + operand.Re;
        res.Im = Im + operand.Im;
        return res;
    }
    complexNumber operator*(complexNumber operand) {
        complexNumber res;
        res.Re = Re * operand.Re - Im * operand.Im;
        res.Im = Re * operand.Im + Im * operand.Re;
        return res;
    }
};

vector<vector<complexNumber>> operator+(vector<vector<complexNumber>> &A, vector<vector<complexNumber>> &B) {
    int size = A.size();
    vector<vector<complexNumber>> res(size, vector<complexNumber>(size, {0.0, 0.0}));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

vector<vector<complexNumber>> operator*(complexNumber coeff, vector<vector<complexNumber>> &matrix) {
    int size = matrix.size();
    vector<vector<complexNumber>> res(size, vector<complexNumber>(size, {0.0, 0.0}));
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            res[i][j] = coeff * matrix[i][j];
        }
    }
    return res;
}

vector<vector<complexNumber>> operator*(vector<vector<complexNumber>> &A, vector<vector<complexNumber>> &B) {
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

void printMatrix(vector<vector<complexNumber>> &res) {
    int n = res.size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << res[i][j].Re << " + i" << res[i][j].Im << "  ";
        }
        cout << "\n";
    }
}

vector<vector<vector<complexNumber>>> computePowers(vector<vector<complexNumber>> &inputMatrix, int limit) {
    int sz = inputMatrix.size();
    vector<vector<vector<complexNumber>>> res(limit + 1, vector<vector<complexNumber>>(sz, vector<complexNumber>(sz, {0.0, 0.0})));
    vector<vector<complexNumber>> temp = inputMatrix;
    vector<vector<complexNumber>> identity(sz, vector<complexNumber>(sz, {0.0, 0.0}));
    for (int i = 0; i < sz; i++) {
        identity[i][i] = {1.0, 0.0};
    }
    res[0] = identity;
    res[1] = inputMatrix;
    for (int i = 2; i <= limit; i++) {
        temp = temp * inputMatrix;
        res[i] = temp;
    }
    return res;
}

vector<vector<complexNumber>> patersonStockmeyerSerial(vector<vector<complexNumber>> &inputMatrix, vector<complexNumber> &coefficients, int polynomialVariable, int polynomialDegree) {
    int inputMatrixSize = inputMatrix.size();
    int degree = coefficients.size() - 1;

    // Compute powers of A till p
    vector<vector<vector<complexNumber>>> powersOfInputMatrix = computePowers(inputMatrix, polynomialVariable);

    // Declare resultant matrix
    vector<vector<complexNumber>> resultantMatrix(inputMatrixSize, vector<complexNumber>(inputMatrixSize, {0.0, 0.0}));

    // Horner's loop
    for (int q = polynomialDegree - 1; q >= 0; q--) {
        vector<vector<complexNumber>> temp(inputMatrixSize, vector<complexNumber>(inputMatrixSize, {0.0, 0.0}));
        for (int j = 0; j < polynomialVariable; j++) {
            complexNumber coeff = {0.0, 0.0};
            int index = polynomialVariable * q + j;
            if (j <= degree) coeff = coefficients[index];
            vector<vector<complexNumber>> adder = coeff * powersOfInputMatrix[j];
            temp = temp + adder;
        }
        vector<vector<complexNumber>> adder = resultantMatrix * powersOfInputMatrix[polynomialVariable];
        resultantMatrix = adder + temp;
    }

    // Return answer
    return resultantMatrix;
}

int main() {
    // Take Input
    int n;
    cin >> n;
    vector<vector<complexNumber>> A(n, vector<complexNumber>(n, {0.0, 0.0}));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    int d;
    cin >> d;
    vector<complexNumber> coeff(d + 1, {0.0, 0.0});
    for (int i = 0; i <= d; i++) {
        double a, b;
        cin >> a >> b;
        coeff[i] = {a, b};
    }
    vector<vector<complexNumber>> res = patersonStockmeyerSerial(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    printMatrix(res);
    return 0;
}