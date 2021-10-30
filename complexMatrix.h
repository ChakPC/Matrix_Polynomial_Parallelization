// Header file for complex matrix operations

#include <bits/stdc++.h>
using namespace std;

inline vector<vector<complex<double>>> operator+(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &B){
    int size = A.size();
    vector<vector<complex<double>>> res(size, vector<complex<double>>(size, complex<double>(0.0, 0.0)));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

inline vector<vector<complex<double>>> operator*(complex<double> coeff, vector<vector<complex<double>>> &matrix){
    int size = matrix.size();
    vector<vector<complex<double>>> res(size, vector<complex<double>>(size, {0.0, 0.0}));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            res[i][j] = coeff * matrix[i][j];
        }
    }
    return res;
}

inline vector<vector<complex<double>>> operator*(vector<vector<complex<double>>> &A, vector<vector<complex<double>>> &B){
    int size = A.size();
    vector<vector<complex<double>>> res(size, vector<complex<double>>(size, {0.0, 0.0}));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                res[i][j] = res[i][j] + A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

inline void printMatrix(vector<vector<complex<double>>> &res){
    int n = res.size();
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            cout << "( " << real(res[i][j]) <<" ) + i( "<< imag(res[i][j]) <<" )  ";
        }
        cout << "\n";
    }
}
