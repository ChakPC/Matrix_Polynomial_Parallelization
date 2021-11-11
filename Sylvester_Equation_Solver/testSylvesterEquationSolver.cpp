#include "../complexMatrix.hpp"
#include "sylvesterEquationSolverSerial.hpp"
#include "sylvesterEquationSolverParallel.hpp"

// Driver function
int main() {
    freopen("../inputs/sylvesterEquationSolver.txt", "r", stdin);
    int p;
    cin >> p;
    complexMatrix A(p, vector<complexNumber>(p, complexNumber(0, 0)));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            double a, b;
            cin >> a >> b;
            A[i][j] = complexNumber(a, b);
        }
    }
    int q;
    cin >> q;
    complexMatrix B(q, vector<complexNumber>(q, complexNumber(0, 0)));
    for (int i = 0; i < q; i++) {
        for (int j = 0; j < q; j++) {
            double a, b;
            cin >> a >> b;
            B[i][j] = complexNumber(a, b);
        }
    }
    complexMatrix C(p, vector<complexNumber>(q, complexNumber(0, 0)));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < q; j++) {
            double a, b;
            cin >> a >> b;
            C[i][j] = complexNumber(a, b);
        }
    }
    complexMatrix X1 = serialSylvesterEquationSolver(A, B, C);
    printMatrix(X1);
    complexMatrix X2 = parallelSylvesterEquationSolver(A, B, C);
    printMatrix(X2);
    return 0;
}