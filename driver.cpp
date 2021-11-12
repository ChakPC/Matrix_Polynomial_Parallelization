#include "complexMatrix.hpp"
#include "complexMatrixSerial.hpp"
#include "complexMatrixParallel.hpp"
#include "Paterson_Stockmeyer/patersonStockmeyerSerial.hpp"
#include "Paterson_Stockmeyer/patersonStockmeyerParallel.hpp"
#include "Schur_Decomposition/schurDecompositionSerial.hpp"
#include "Schur_Decomposition/schurDecompositionParallel.hpp"
#include "Sylvester_Equation_Solver/sylvesterEquationSolverSerial.hpp"
#include "Sylvester_Equation_Solver/sylvesterEquationSolverParallel.hpp"
#include "Parlett_Recurrence/parletRecurrenceParallel.hpp"
#include "Parlett_Recurrence/parletRecurrenceSerial.hpp"

#define matrixSizeLimit 5
#define polynomialDegreeLimit 5

int main() {
    freopen("./Analytics/patersonStockmeyer/input.txt", "r", stdin);
    // Take Input
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

    if (A.size() <= matrixSizeLimit || d <= polynomialDegreeLimit) {
        // Output from Paterson Stockmeyer Serial
        complexMatrix patersonStockmeyerSerialResult = patersonStockmeyerSerial(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
        printMatrix(patersonStockmeyerSerialResult);

        // Output from Paterson Stockmeyer Parallel
        complexMatrix patersonStockmeyerParallelResult = patersonStockmeyerParallel(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
        printMatrix(patersonStockmeyerParallelResult);
    } else {
        // Output from Parlett Recurrence Serial
        complexMatrix parlettRecurrenceSerialResult = parlettRecurrenceSerial(A, coeff, 10, sqrt(d) + 1, sqrt(d) + 1);
        printMatrix(parlettRecurrenceSerialResult);

        // Output from Parlett Recurrence Parallel
        complexMatrix parlettRecurrenceParallelResult = parlettRecurrenceParallel(A, coeff, 10, sqrt(d) + 1, sqrt(d) + 1);
        printMatrix(parlettRecurrenceParallelResult);
    }

    return 0;
}