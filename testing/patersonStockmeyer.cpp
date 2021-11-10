#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "../complexMatrix.hpp"
#include "../Paterson_Stockmeyer/patersonStockmeyer.hpp"

TEST_CASE("Paterson Stockmeyer Serial", "") {
    complexMatrix A = {
        {complexNumber(1, 0), complexNumber(0, 0)},
        {complexNumber(0, 0), complexNumber(2, 0)}};

    vector<complexNumber> coeff = {
        complexNumber(1, 0),
        complexNumber(2, 0),
        complexNumber(3, 0)};

    complexMatrix res = {
        {complexNumber(6, 0), complexNumber(0, 0)},
        {complexNumber(0, 0), complexNumber(17, 0)}};

    REQUIRE(res == patersonStockmeyerSerial(A, coeff, 2, 2));
}

TEST_CASE("Paterson Stockmeyer Parallel", "") {
    complexMatrix A = {
        {complexNumber(1, 0), complexNumber(0, 0)},
        {complexNumber(0, 0), complexNumber(2, 0)}};

    vector<complexNumber> coeff = {
        complexNumber(1, 0),
        complexNumber(2, 0),
        complexNumber(3, 0)};

    complexMatrix res = {
        {complexNumber(6, 0), complexNumber(0, 0)},
        {complexNumber(0, 0), complexNumber(17, 0)}};

    REQUIRE(res == patersonStockmeyerParallel(A, coeff, 2, 2));
}