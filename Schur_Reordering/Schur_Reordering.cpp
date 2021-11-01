#include "../complexMatrix.h"
#define ZERO 1e-10

bool realCmp(pair<complexNumber, int> &a, pair<complexNumber, int> &b) {
    double z1 = real(a.first) * real(a.first) + imag(a.first) * imag(a.first);
    double z2 = real(b.first) * real(b.first) + imag(b.first) * imag(b.first);
    if (imag(a.first) == 0 && imag(b.first) == 0) {
        return real(a.first) < real(b.first);
    }
    return z1 < z2;
}

bool imagCmp(pair<complexNumber, int> &a, pair<complexNumber, int> &b) {
}

void schurReorderingSerial(complexMatrix &inputMatrix) {
    int n = inputMatrix.size();
    vector<pair<complexNumber, int>> realEigenValues;
    vector<pair<vector<complexNumber>, int>> imaginaryEigenValues;
    for (int i = 0; i < n; i++) {
        if (i < n - 1 && abs(inputMatrix[i + 1][i]) > ZERO) {
            imaginaryEigenValues.push_back({{inputMatrix[i][i], inputMatrix[i][i + 1], inputMatrix[i + 1][i], inputMatrix[i + 1][i + 1]}, i});
            i++;
        } else {
            realEigenValues.push_back({inputMatrix[i][i], i});
        }
    }
    sort(realEigenValues.begin(), realEigenValues.end(), realCmp);
    sort(imaginaryEigenValues.begin(), imaginaryEigenValues.end(), imagCmp);
}
int main() {
    return 0;
}
