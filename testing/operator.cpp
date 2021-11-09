#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../complexMatrix.hpp"

TEST_CASE("Operator Overloading", "") {
    complexMatrix A(2, vector<complexNumber>(2, complexNumber(0, 0)));
    complexMatrix B(2, vector<complexNumber>(2, complexNumber(0, 0)));

    A = {
        {complexNumber(1, 3), complexNumber(2, 5)},
        {complexNumber(8, 0), complexNumber(10, 7)},
    };
    B = {
        {complexNumber(55, 2), complexNumber(12, 2)},
        {complexNumber(1, 3), complexNumber(34, 3)},
    };

    complexMatrix sum = {
        {complexNumber(56, 5), complexNumber(14, 7)},
        {complexNumber(9, 3), complexNumber(44, 10)}};

    complexMatrix product = {
        {complexNumber(36, 178), complexNumber(59, 214)},
        {complexNumber(429, 53), complexNumber(415, 284)}};

    complexMatrix twoB = {
        {complexNumber(110, 4), complexNumber(24, 4)},
        {complexNumber(2, 6), complexNumber(68, 6)}};

    REQUIRE(sum == A + B);
    REQUIRE(product == A * B);
    REQUIRE(twoB == 2 * B);
}
