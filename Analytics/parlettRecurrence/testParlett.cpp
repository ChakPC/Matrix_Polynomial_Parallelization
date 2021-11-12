#include "../../complexMatrix.hpp"
#include "../../complexMatrixSerial.hpp"
#include "../../complexMatrixParallel.hpp"
#include "../../Paterson_Stockmeyer/patersonStockmeyerSerial.hpp"
#include "../../Paterson_Stockmeyer/patersonStockmeyerParallel.hpp"
#include "../../Schur_Decomposition/schurDecompositionSerial.hpp"
#include "../../Schur_Decomposition/schurDecompositionParallel.hpp"
#include "../../Sylvester_Equation_Solver/sylvesterEquationSolverSerial.hpp"
#include "../../Sylvester_Equation_Solver/sylvesterEquationSolverParallel.hpp"
#include "../../Parlett_Recurrence/parletRecurrenceParallel.hpp"
#include "../../Parlett_Recurrence/parletRecurrenceSerial.hpp"

// Driver function
int main() {
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "a", stdout);

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
    complexMatrix parlettRecurrenceSerialResult = parlettRecurrenceSerial(A, coeff, 10, sqrt(d) + 1, sqrt(d) + 1);
    printf("S: %f\n", omp_get_wtime() - start);

    omp_set_num_threads(8);
    start = omp_get_wtime();
    complexMatrix parlettRecurrenceParallelResult = parlettRecurrenceParallel(A, coeff, 10, sqrt(d) + 1, sqrt(d) + 1);
    printf("P: %f\n", omp_get_wtime() - start);

    return 0;
}