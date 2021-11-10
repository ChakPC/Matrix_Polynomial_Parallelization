#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../complexMatrix.hpp"

std::string FullPrecision(double d) {
    auto s = std::ostringstream{};
    s << std::setprecision(std::numeric_limits<double>::max_digits10) << d;
    return s.str();
}

TEST_CASE("Operator Overloading", "") {
    complexMatrix A = {
        {complexNumber(1, 3), complexNumber(2, 5)},
        {complexNumber(8, 0), complexNumber(10, 7)},
    };
    complexMatrix A1 = {
        {complexNumber(1.000001, 3.000000013), complexNumber(2, 5)},
        {complexNumber(8, 0), complexNumber(10, 7)},
    };
    complexMatrix B = {
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

    REQUIRE(areMatricesEqual(A, A1) == true);
    REQUIRE(sum == A + B);
    REQUIRE(product == A * B);
    REQUIRE(twoB == 2 * B);
}
