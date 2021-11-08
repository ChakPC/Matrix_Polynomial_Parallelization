// Header file for serial Implementation of Sylvester Equation Solver

// Function to compute K = Iq * A - B.T * Ip
// * denotes kronecker product
inline complexMatrix computeK(complexMatrix &A, complexMatrix &B) {
    int p = A.size();
    int q = B.size();
    complexMatrix Iq = identityMatrix(q);
    complexMatrix Ip = identityMatrix(p);
    complexMatrix t1 = computeKroneckerProduct(Iq, A);
    complexMatrix BTranspose = transposeMatrix(B);
    BTranspose = complexNumber(-1, 0) * BTranspose;
    complexMatrix t2 = computeKroneckerProduct(BTranspose, Ip);
    complexMatrix res = t1 + t2;
    return res;
}

// Function to solve sylvester equation AX - XB = C
// A: sqaure matrix (p x p)
// B: square matrix (q x q)
// X: (p x q) matrix
// C: (p x q) matrix
inline complexMatrix sylvesterEquationSolver(complexMatrix A, complexMatrix B, complexMatrix C) {
    int p = A.size(), q = B.size();
    complexMatrix x(p * q, vector<complexNumber>(1, complexNumber(0, 0)));
    complexMatrix K = computeK(A, B);
    complexMatrix b = columnStack(C);
    x = gaussElimination(K, b);
    complexMatrix X = reverseStacking(x, p);
    processZero(X);
    return X;
}

