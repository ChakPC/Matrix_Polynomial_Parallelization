#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../complexMatrix.hpp"
#include "../complexMatrixSerial.hpp"
#include "../complexMatrixParallel.hpp"
#include "../Schur_Decomposition/schurDecompositionSerial.hpp"
#include "../Schur_Decomposition/schurDecompositionParallel.hpp"

TEST_CASE("Schur Decomposition Random Test", "") {
    system("python3 testGeneratorSchurDecomposition.py 2");
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

    omp_set_num_threads(1);
    double start = omp_get_wtime();
    auto resSerial = schurDecompositionSerial(A, 100);
    cout << "Serial:\n";
    printf("S: %f\n", omp_get_wtime() - start);
    printMatrix(resSerial[0]);
    printMatrix(resSerial[1]);

    start = omp_get_wtime();
    omp_set_num_threads(8);
    auto resParallel = schurDecompositionParallel(A, 100);
    cout << "Paralel:\n";
    printf("P: %f\n", omp_get_wtime() - start);
    printMatrix(resParallel[0]);
    printMatrix(resParallel[1]);
    // REQUIRE(areMatricesEqual(resSerial, resParallel) == true);
}