#include <bits/stdc++.h>
// #include "/usr/lib/gcc/x86_64-linux-gnu/9/include/omp.h"
#include <omp.h>

using namespace std;
#define complexNumber complex<double>
#define complexMatrix vector<vector<complexNumber>>
#define clusterThreshold 0.1
#define zeroLimit 1e-5

inline bool areSame(complexNumber x, complexNumber y) {
    string realX = to_string((long long)real(x));
    string imagX = to_string((long long)imag(x));

    string realY = to_string((long long)real(y));
    string imagY = to_string((long long)imag(y));

    if (realX.size() != realY.size() || imagX.size() != imagY.size()) return false;

    bool flag = true;
    int sz = min(3, (int)realX.size());
    for (int i = 0; i < sz; i++) {
        if (realX[i] != realY[i] || imagX[i] != imagY[i]) flag = false;
    }
    if (!flag) return false;
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

// Comparator to compare complex number
inline bool cmp(complexNumber &A, complexNumber &B) {
    double a = real(A) * real(A) + imag(A) * imag(A);
    double b = real(B) * real(B) + imag(B) * imag(B);
    return a < b;
}

/*
 * check if two matrices are equal by checking if the difference between two
 * values is smaller than zeroLimit
 * @param A Matrix A
 * @param B Matrix B
 * @return true if A == B else false
 */
inline bool areMatricesEqual(complexMatrix &A, complexMatrix &B) {
    if (A.size() != B.size() || A[0].size() != B[0].size()) return false;
    volatile bool flag = false;
    // complexNumber x, y;
#pragma omp parallel for shared(flag) collapse(2)
    for (int i = 0; i < A.size(); i++) {
        for (int j = 0; j < A[0].size(); j++) {
            if (flag) continue;
            if (!areSame(A[i][j], B[i][j])) {
                // x = A[i][j];
                // y = B[i][j];
                flag = true;
            }
        }
    }
    // printComplex(x);
    // printComplex(y);
    if (flag) return false;
    return true;
}
