// Serial Implementation of Paterson - Stockmeyer algorithm
#include "../complexMatrix.h"
using namespace std;

// Function to compute powers of input matrix till exponent limit
// Brute Force Implmentation: O(limit * n * n)
vector<complexMatrix> computePowers(complexMatrix &inputMatrix, int limit){
    int sz = inputMatrix.size();
    vector<complexMatrix> res(limit+1, complexMatrix(sz, vector<complexNumber>(sz, complexNumber(0.0, 0.0))));
    complexMatrix temp = inputMatrix;
    complexMatrix identity(sz, vector<complexNumber>(sz, complexNumber(0.0, 0.0)));
    for(int i = 0; i < sz; i++){
        identity[i][i] = complexNumber(1.0, 0.0);
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
complexMatrix patersonStockmeyerSerial(complexMatrix &inputMatrix, vector<complexNumber> &coefficients, int polynomialVariable, int polynomialDegree){
    int inputMatrixSize = inputMatrix.size();
    int degree = coefficients.size() - 1;

    // Compute powers of A till p
    vector<complexMatrix> powersOfInputMatrix = computePowers(inputMatrix, polynomialVariable);

    // Declare resultant matrix
    complexMatrix resultantMatrix(inputMatrixSize, vector<complexNumber>(inputMatrixSize, complexNumber(0.0, 0.0)));

    // Horner's loop
    for(int q = polynomialDegree-1; q >= 0; q--){
        complexMatrix temp(inputMatrixSize, vector<complexNumber>(inputMatrixSize, complexNumber(0.0, 0.0)));
        for(int j = 0; j < polynomialVariable; j++){
            complexNumber coeff = {0.0, 0.0};
            int index = polynomialVariable * q + j;
            if(j <= degree) coeff = coefficients[index];
            complexMatrix adder = coeff * powersOfInputMatrix[j];
            temp = temp + adder;
        }
        complexMatrix adder = resultantMatrix * powersOfInputMatrix[polynomialVariable];
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
    complexMatrix A(n, vector<complexNumber>(n, {0.0, 0.0}));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    int d;
    cin >> d;
    vector<complexNumber> coeff(d+1, {0.0, 0.0});
    for(int i=0; i<=d; i++){
        double a, b;
        cin >> a >> b;
        coeff[i] = {a, b};
    }
    complexMatrix res = patersonStockmeyerSerial(A, coeff, sqrt(d)+1, sqrt(d)+1);
    printMatrix(res);
    return 0;
}