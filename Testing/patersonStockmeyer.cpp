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

    complexMatrix output = patersonStockmeyerSerial(A, coeff, 2, 2);
    REQUIRE(areMatricesEqual(res, output) == true);
}

TEST_CASE("Paterson Stockmeyer Parallel", "Checking with 2d Matrices") {
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

    complexMatrix patersonStockmeyerParallelOutput = patersonStockmeyerParallel(A, coeff, 2, 2);
    REQUIRE(areMatricesEqual(res, patersonStockmeyerParallelOutput) == true);
}

TEST_CASE("Paterson Stockmeyer Parallel (complex)", "") {
    complexMatrix A = {
        {complexNumber(2, 0), complexNumber(1, 1), complexNumber(4, 0)},
        {complexNumber(0, 0), complexNumber(5, 0), complexNumber(0, 4)},
        {complexNumber(8, 0), complexNumber(1, 0), complexNumber(0, 5)}};

    vector<complexNumber> coeff = {
        complexNumber(4, 0),
        complexNumber(0, 5),
        complexNumber(10, 0)};

    complexMatrix res = {
        {complexNumber(364, 10), complexNumber(105, 75), complexNumber(40, 260)},
        {complexNumber(0, 320), complexNumber(254, 65), complexNumber(-220, 200)},
        {complexNumber(160, 440), complexNumber(130, 135), complexNumber(49, 40)}};

    complexMatrix patersonStockmeyerParallelOutput3d = patersonStockmeyerParallel(A, coeff, 2, 2);
    REQUIRE(areMatricesEqual(res, patersonStockmeyerParallelOutput3d) == true);
}