#include "../complexMatrix.hpp"

int main() {
    complexMatrix A = {
        {complexNumber(1, 3), complexNumber(2, 5)},
        {complexNumber(8, 0), complexNumber(10, 7)},
    };
    complexMatrix B = {
        {complexNumber(1.0001, 3.000000013), complexNumber(2, 5)},
        {complexNumber(8, 0), complexNumber(10, 7)},
    };

    if (A == B) {
        cout << "Equal";
    } else {
        cout << "Not Equal";
    }
}