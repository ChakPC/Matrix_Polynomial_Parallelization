#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../complexMatrix.hpp"
#include "../complexMatrixParallel.hpp"
#include "../complexMatrixSerial.hpp"
#include "../Schur_Decomposition/schurDecompositionSerial.hpp"
#include "../Schur_Decomposition/schurDecompositionParallel.hpp"
#include "../Sylvester_Equation_Solver/sylvesterEquationSolverSerial.hpp"
#include "../Sylvester_Equation_Solver/sylvesterEquationSolverParallel.hpp"
#include "../Paterson_Stockmeyer/patersonStockmeyerParallel.hpp"
#include "../Paterson_Stockmeyer/patersonStockmeyerSerial.hpp"
#include "../Parlett_Recurrence/parletRecurrenceParallel.hpp"
#include "../Parlett_Recurrence/parletRecurrenceSerial.hpp"

TEST_CASE("Paterson Stockmeyer stress test", "") {
    system("python3 testGeneratorPatersonStockmeyer.py 50");
    freopen("input.txt", "r", stdin);
    int n;
    cin >> n;
    complexMatrix A(n, vector<complexNumber>(n, {0.0, 0.0}));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    int d;
    cin >> d;
    vector<complexNumber> coeff(d + 1, {0.0, 0.0});
    for (int i = 0; i <= d; i++) {
        double a, b;
        cin >> a >> b;
        coeff[i] = {a, b};
    }
    omp_set_num_threads(1);
    double start = omp_get_wtime();
    complexMatrix resSerial = parlettRecurrenceSerial(A, coeff, 10, sqrt(d) + 1, sqrt(d) + 1);
    printf("S: %f\n", omp_get_wtime() - start);
    start = omp_get_wtime();

    omp_set_num_threads(8);
    complexMatrix resParallel = parlettRecurrenceParallel(A, coeff, 10, sqrt(d) + 1, sqrt(d) + 1);
    printf("P: %f\n", omp_get_wtime() - start);
    REQUIRE(areMatricesEqual(resSerial, resParallel) == true);
}