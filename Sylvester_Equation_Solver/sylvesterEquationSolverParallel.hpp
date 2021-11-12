/*
 * Function to compute K
 * denotes kronecker product
 * @param A Matrix A
 * @param B Matrix B
 * @return K (Iq * A - B.T * Ip, where * is Kronecker product)
 */
inline complexMatrix computeKParallel(complexMatrix &A, complexMatrix &B) {
    int p = A.size();
    int q = B.size();
    complexMatrix Iq = identityMatrixParallel(q);
    complexMatrix Ip = identityMatrixParallel(p);
    complexMatrix t1 = computeKroneckerProductParallel(Iq, A);
    complexMatrix BTranspose = transposeMatrixParallel(B);
    BTranspose = complexNumber(-1, 0) * BTranspose;
    complexMatrix t2 = computeKroneckerProductParallel(BTranspose, Ip);
    complexMatrix res = t1 + t2;
    return res;
}

/*
 * Function to solve sylvester equation AX - XB = C
 * @param A square Matrix (p x p)
 * @param B square Matrix (q x q)
 * @param C (p x q) matrix
 * @return X, result of sylvester Equation
 */
inline complexMatrix sylvesterEquationSolverParallel(complexMatrix A, complexMatrix B, complexMatrix C) {
    int p = A.size(), q = B.size();
    complexMatrix x(p * q, vector<complexNumber>(1, complexNumber(0, 0)));
    complexMatrix K = computeKParallel(A, B);
    complexMatrix b = columnStackParallel(C);
    x = gaussEliminationParallel(K, b);
    complexMatrix X = reverseStackingParallel(x, p);
    processZeroParallel(X);
    return X;
}
