#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../complexMatrix.hpp"
#include "../complexMatrixSerial.hpp"
#include "../complexMatrixParallel.hpp"
#include "../Sylvester_Equation_Solver/sylvesterEquationSolverSerial.hpp"
#include "../Sylvester_Equation_Solver/sylvesterEquationSolverParallel.hpp"

TEST_CASE("Sylvester Equation Solver", "") {
    system("python3 testSylvesterEquationSolver.py");
    freopen("input.txt", "r", stdin);
    int p;
    cin >> p;
    complexMatrix A(p, vector<complexNumber>(p, {0.0, 0.0}));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    int q;
    cin >> q;
    complexMatrix B(p, vector<complexNumber>(p, {0.0, 0.0}));
    for (int i = 0; i < p; i++) {
        for (int j = 0; j < p; j++) {
            double a, b;
            cin >> a >> b;
            B[i][j] = {a, b};
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

    omp_set_num_threads(1);
    double start = omp_get_wtime();
    complexMatrix resSerial = sylvesterEquationSolverSerial(A, B, C);
    printf("S: %f\n", omp_get_wtime() - start);

    omp_set_num_threads(8);
    start = omp_get_wtime();
    complexMatrix resParallel = sylvesterEquationSolverParallel(A, B, C);
    printf("P: %f\n", omp_get_wtime() - start);

    REQUIRE(areMatricesEqual(resSerial, resParallel) == true);
}