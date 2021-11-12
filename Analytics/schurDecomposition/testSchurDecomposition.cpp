#include "../../complexMatrix.hpp"
#include "../../complexMatrixSerial.hpp"
#include "../../complexMatrixParallel.hpp"
#include "../../Schur_Decomposition/schurDecompositionSerial.hpp"
#include "../../Schur_Decomposition/schurDecompositionParallel.hpp"

// Driver function for testing
int main() {
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "a", stdout);
    // Take Input
    int n;
    cin >> n;
    complexMatrix A(n, vector<complexNumber>(n, {0.0, 0.0}));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a, b;
            cin >> a >> b;
            A[i][j] = complexNumber(a, b);
        }
    }

    omp_set_num_threads(1);
    double start = omp_get_wtime();
    vector<complexMatrix> t1 = schurDecompositionSerial(A, 10);
    printf("S: %f\n", omp_get_wtime() - start);

    omp_set_num_threads(8);
    start = omp_get_wtime();
    vector<complexMatrix> t2 = schurDecompositionParallel(A, 10);
    printf("P: %f\n", omp_get_wtime() - start);
    return 0;
}