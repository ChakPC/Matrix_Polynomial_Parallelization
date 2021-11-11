#include "../complexMatrix.hpp"
#include "patersonStockmeyerSerial.hpp"
#include "patersonStockmeyerParallel.hpp"

// Driver function for testing
int main() {
    // Take Input
    freopen("../inputs/patersonStockmeyer.txt", "r", stdin);
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
    complexMatrix res1 = serialPatersonStockmeyer(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    printMatrix(res1);
    complexMatrix res2 = parallelPatersonStockmeyer(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    printMatrix(res2);
    return 0;
}
