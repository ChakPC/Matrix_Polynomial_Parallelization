#include "../../complexMatrix.hpp"
#include "../../complexMatrixSerial.hpp"
#include "../../complexMatrixParallel.hpp"
#include "../../Paterson_Stockmeyer/patersonStockmeyerSerial.hpp"
#include "../../Paterson_Stockmeyer/patersonStockmeyerParallel.hpp"

// Driver function for testing
int main() {
    // Take Input
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
    complexMatrix resSerial = patersonStockmeyerSerial(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    printf("S: %f\n", omp_get_wtime() - start);

    omp_set_num_threads(8);
    start = omp_get_wtime();
    complexMatrix resParallel = patersonStockmeyerParallel(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    printf("P: %f\n", omp_get_wtime() - start);
    return 0;
}
