#include "../complexMatrix.hpp"
#include "../complexMatrixSerial.hpp"
#include "../complexMatrixParallel.hpp"
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
    auto t1 = std::chrono::high_resolution_clock::now();
    complexMatrix res1 = patersonStockmeyerSerial(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    auto t2 = std::chrono::high_resolution_clock::now();
    printMatrix(res1);
    auto t3 = std::chrono::high_resolution_clock::now();
    complexMatrix res2 = patersonStockmeyerParallel(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    auto t4 = std::chrono::high_resolution_clock::now();
    printMatrix(res2);
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
    cout << "Serial Time: " << duration.count() << endl;
    duration = std::chrono::duration_cast<std::chrono::microseconds>(t4 - t3);
    cout << "Parallel Time: " << duration.count() << endl;
    return 0;
}
