// Serial Implementation of Paterson - Stockmeyer algorithm
#include "../complexMatrix.h"
using namespace std;

// Function to compute powers of input matrix till exponent limit
// Brute Force Implmentation: O(limit * n * n)
vector<vector<vector<complex<double>>>> computePowers(vector<vector<complex<double>>> &inputMatrix, int limit){
    int sz = inputMatrix.size();
    vector<vector<vector<complex<double>>>> res(limit+1, vector<vector<complex<double>>>(sz, vector<complex<double>>(sz, complex<double>(0.0, 0.0))));
    vector<vector<complex<double>>> temp = inputMatrix;
    vector<vector<complex<double>>> identity(sz, vector<complex<double>>(sz, complex<double>(0.0, 0.0)));
    for(int i = 0; i < sz; i++){
        identity[i][i] = complex<double>(1.0, 0.0);
    }
    res[0] = identity;
    res[1] = inputMatrix;
    for (int i = 2; i <= limit; i++) {
        temp = temp * inputMatrix;
        res[i] = temp;
    }
    return res;
}

// Serial Implementation of Paterson Stockmeyer
vector<vector<complex<double>>> patersonStockmeyerSerial(vector<vector<complex<double>>> &inputMatrix, vector<complex<double>> &coefficients, int polynomialVariable, int polynomialDegree){
    int inputMatrixSize = inputMatrix.size();
    int degree = coefficients.size() - 1;

    // Compute powers of A till p
    vector<vector<vector<complex<double>>>> powersOfInputMatrix = computePowers(inputMatrix, polynomialVariable);

    // Declare resultant matrix
    vector<vector<complex<double>>> resultantMatrix(inputMatrixSize, vector<complex<double>>(inputMatrixSize, complex<double>(0.0, 0.0)));

    // Horner's loop
    for(int q = polynomialDegree-1; q >= 0; q--){
        vector<vector<complex<double>>> temp(inputMatrixSize, vector<complex<double>>(inputMatrixSize, complex<double>(0.0, 0.0)));
        for(int j = 0; j < polynomialVariable; j++){
            complex<double> coeff = {0.0, 0.0};
            int index = polynomialVariable * q + j;
            if(j <= degree) coeff = coefficients[index];
            vector<vector<complex<double>>> adder = coeff * powersOfInputMatrix[j];
            temp = temp + adder;
        }
        vector<vector<complex<double>>> adder = resultantMatrix * powersOfInputMatrix[polynomialVariable];
        resultantMatrix = adder + temp;
    }

    // Return answer
    return resultantMatrix;
}

// Driver function for testing
int main() {
    // Take Input
    int n;
    cin >> n;
    vector<vector<complex<double>>> A(n, vector<complex<double>>(n, {0.0, 0.0}));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    int d;
    cin >> d;
    vector<complex<double>> coeff(d+1, {0.0, 0.0});
    for(int i=0; i<=d; i++){
        double a, b;
        cin >> a >> b;
        coeff[i] = {a, b};
    }
    vector<vector<complex<double>>> res = patersonStockmeyerSerial(A, coeff, sqrt(d)+1, sqrt(d)+1);
    printMatrix(res);
    return 0;
}