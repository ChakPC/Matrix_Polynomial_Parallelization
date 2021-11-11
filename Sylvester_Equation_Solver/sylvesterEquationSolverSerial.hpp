// Header file for serial Implementation of Sylvester Equation Solver

// Function to compute K = Iq * A - B.T * Ip
// * denotes kronecker product
inline complexMatrix computeKSerial(complexMatrix &A, complexMatrix &B) {
    int p = A.size();
    int q = B.size();
    complexMatrix Iq = identityMatrixSerial(q);
    complexMatrix Ip = identityMatrixSerial(p);
    complexMatrix t1 = computeKroneckerProductSerial(Iq, A);
    complexMatrix BTranspose = transposeMatrixSerial(B);
    BTranspose = multiplySerial(complexNumber(-1, 0), BTranspose);
    complexMatrix t2 = computeKroneckerProductSerial(BTranspose, Ip);
    complexMatrix res = addSerial(t1, t2);
    return res;
}

// Function to solve sylvester equation AX - XB = C
// A: sqaure matrix (p x p)
// B: square matrix (q x q)
// X: (p x q) matrix
// C: (p x q) matrix
inline complexMatrix sylvesterEquationSolverSerial(complexMatrix A, complexMatrix B, complexMatrix C) {
    int p = A.size(), q = B.size();
    complexMatrix x(p * q, vector<complexNumber>(1, complexNumber(0, 0)));
    complexMatrix K = computeKSerial(A, B);
    complexMatrix b = columnStackSerial(C);
    x = gaussEliminationSerial(K, b);
    complexMatrix X = reverseStackingSerial(x, p);
    processZeroSerial(X);
    return X;
}
