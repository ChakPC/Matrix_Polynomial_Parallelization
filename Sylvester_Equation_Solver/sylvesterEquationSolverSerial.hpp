/*
 * Function to compute K =
 * denotes kronecker product
 * @param A Matrix A
 * @param B Matrix B
 * @return K (Iq * A - B.T * Ip, where * is Kronecker product)
 */
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

/*
 * Function to solve sylvester equation AX - XB = C
 * @param A square Matrix (p x p)
 * @param B square Matrix (q x q)
 * @param C (p x q) matrix
 * @return X, result of sylvester Equation
 */
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
