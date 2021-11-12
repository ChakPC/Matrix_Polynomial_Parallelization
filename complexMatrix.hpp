#include <bits/stdc++.h>
#include "/usr/lib/gcc/x86_64-linux-gnu/9/include/omp.h"
// #include <omp.h>

using namespace std;
#define complexNumber complex<double>
#define complexMatrix vector<vector<complexNumber>>
#define clusterThreshold 0.1
#define zeroLimit 1e-5

inline bool areSame(complexNumber x, complexNumber y) {
    return (fabs(real(x) - real(y)) < zeroLimit) && (fabs(imag(x) - imag(y)) < zeroLimit);
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

// Comparator to compare complex number
inline bool cmp(complexNumber &A, complexNumber &B) {
    double a = real(A) * real(A) + imag(A) * imag(A);
    double b = real(B) * real(B) + imag(B) * imag(B);
    return a < b;
}