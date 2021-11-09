#include "../complexMatrix.hpp"
#include "schurDecompositionSerial.hpp"
#include <omp.h>
// Driver function for testing
int main() {
    omp_set_num_threads(4);
    freopen("../inputs/schurDecomposition.txt", "r", stdin);
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
    vector<complexMatrix> t = schurDecompositionSerial(A, 10);
    printf("\n Q: \n");
    printMatrix(t[0]);
    printf("\n Ak: \n");
    printMatrix(t[1]);
    return 0;
}