#include "complexMatrix.hpp"
#include "Paterson_Stockmeyer/patersonStockmeyerSerial.hpp"
#include "Schur_Decomposition/schurDecompositionSerial.hpp"
#include "Sylvester_Equation_Solver/sylvesterEquationSolverSerial.hpp"

// Function to construct blocks in quasi upper triangular matrix based on diagonal values
tuple<int, map<pair<int, int>, complexMatrix>, map<pair<int, int>, pair<int, int>>> constructBlocks(complexMatrix &inputMatrix) {
    // res: to store block index and block matrix
    // limits: to store starting row and starting column of a block of given block index
    map<pair<int, int>, complexMatrix> res;
    map<pair<int, int>, pair<int, int>> limits;

    // Store diagonal blocks
    int size = inputMatrix.size();
    int count = 0;
    for (int i = 0; i < size; i++) {
        count++;
        if (i + 1 < size && ((abs(real(inputMatrix[i + 1][i])) > zeroLimit) || (abs(imag(inputMatrix[i + 1][i])) > zeroLimit))) {
            complexMatrix block(2, vector<complexNumber>(2));
            block[0][0] = inputMatrix[i][i];
            block[0][1] = inputMatrix[i][i + 1];
            block[1][0] = inputMatrix[i + 1][i];
            block[1][1] = inputMatrix[i + 1][i + 1];
            res[{count, count}] = block;
            limits[{count, count}] = {i, i};
            i++;
        } else {
            complexMatrix block(1, vector<complexNumber>(1));
            block[0][0] = inputMatrix[i][i];
            res[{count, count}] = block;
            limits[{count, count}] = {i, i};
        }
    }

    // Store off-diagonal blocks
    for (int i = 1; i <= count; i++) {
        for (int j = i + 1; j <= count; j++) {
            int row = res[{i, i}].size();
            int col = res[{j, j}].size();
            complexMatrix block(row, vector<complexNumber>(col));
            int rowStart = limits[{i, i}].first;
            int colStart = limits[{j, j}].second;
            for (int x = rowStart; x < rowStart + row; x++) {
                for (int y = colStart; y < colStart + col; y++) {
                    block[x - rowStart][y - colStart] = inputMatrix[x][y];
                }
            }
            res[{i, j}] = block;
        }
    }

    // count: No. of diagonal blocks
    // res: map of block index and corresponding block
    return {count, res, limits};
}

// Function to compute matrix polynomial using Parlett recurrence
complexMatrix parlettRecurrenceSerial(complexMatrix &A, vector<complexNumber> &coeff, int numOfIterations, int PS_p, int PS_s) {
    // Schur Decomposition
    vector<complexMatrix> schur = schurDecompositionSerial(A, numOfIterations);
    if(schur.size() == 0){
        // Check compatibility for Given's rotation
        complexMatrix patersonStockmeyerResult = patersonStockmeyerSerial(A, coeff, PS_p, PS_s);
        return patersonStockmeyerResult;
    }
    complexMatrix Q = schur[0];
    complexMatrix T = schur[1];

    // Check clustering of Eigen values
    bool flag = checkClusteredEigenValues(T);
    if(flag){
        complexMatrix patersonStockmeyerResult = patersonStockmeyerSerial(A, coeff, PS_p, PS_s);
        return patersonStockmeyerResult;
    }

    // Divide matri xinto blocks
    map<pair<int, int>, complexMatrix> computedBlocks;
    auto res = constructBlocks(T);
    int numOfBlocks = get<0>(res);
    map<pair<int, int>, complexMatrix> blocks = get<1>(res);
    map<pair<int, int>, pair<int, int>> limits = get<2>(res);

    // Compute function for diagonal blocks using Paterson Stockmeyer
    int size = A.size();
    for (int i = 1; i <= numOfBlocks; i++) {
        complexMatrix computedBlock = patersonStockmeyerSerial(blocks[{i, i}], coeff, PS_p, PS_s);
        processZero(computedBlock);
        computedBlocks[{i, i}] = computedBlock;
    }

    // Compute function for off-diagonal blocks in diagonal fashion using Sylvester Equation Solver
    for (int diff = 1; diff < numOfBlocks; diff++) {
        for (int i = 1; i + diff <= numOfBlocks; i++) {
            int j = i + diff;
            complexMatrix SES_A = blocks[{i, i}];
            complexMatrix SES_B = blocks[{j, j}];
            complexMatrix SES_C = computedBlocks[{i, i}] * blocks[{i, j}];
            complexMatrix temp = blocks[{i, j}] * computedBlocks[{j, j}];
            temp = complexNumber(-1, 0) * temp;
            SES_C = SES_C + temp;
            for (int k = i + 1; k < j; k++) {
                complexMatrix term1 = computedBlocks[{i, k}] * blocks[{k, j}];
                complexMatrix term2 = blocks[{i, k}] * computedBlocks[{k, j}];
                term2 = complexNumber(-1, 0) * term2;
                term1 = term1 + term2;
                SES_C = SES_C + term1;
            }
            complexMatrix SES_X = sylvesterEquationSolver(SES_A, SES_B, SES_C);
            computedBlocks[{i, j}] = SES_X;
        }
    }

    // Resultant function matrix for quasi upper triangular matrix
    complexMatrix functionOnT(size, vector<complexNumber>(size, complexNumber(0, 0)));
    for (int i = 1; i <= numOfBlocks; i++) {
        for (int j = i; j <= numOfBlocks; j++) {
            int rowStart = limits[{i, i}].first;
            int colStart = limits[{j, j}].second;
            int row = computedBlocks[{i, j}].size();
            int col = computedBlocks[{i, j}][0].size();
            for (int x = rowStart; x < rowStart + row; x++) {
                for (int y = colStart; y < colStart + col; y++) {
                    functionOnT[x][y] = computedBlocks[{i, j}][x - rowStart][y - colStart];
                    if (abs(real(functionOnT[x][y])) < 1e-5 && abs(imag(functionOnT[x][y])) < 1e-5) {
                        functionOnT[x][y] = complexNumber(0, 0);
                    }
                }
            }
        }
    }

    // Resultant function matrix for input matrix
    complexMatrix finalMatrix = Q * functionOnT;
    complexMatrix QTranspose = transposeMatrix(Q);
    finalMatrix = finalMatrix * QTranspose;

    return finalMatrix;
}

int main() {
    freopen("inputs/driver[33].txt", "r", stdin);
    // Take Input
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

    // Output from Paterson Stockmeyer
    complexMatrix patersonStockmeyerResult = patersonStockmeyerSerial(A, coeff, sqrt(d) + 1, sqrt(d) + 1);
    printMatrix(patersonStockmeyerResult);

    // Output from Parlett Recurrence
    complexMatrix parlettRecurrenceResult = parlettRecurrenceSerial(A, coeff, 10, sqrt(d) + 1, sqrt(d) + 1);
    printMatrix(parlettRecurrenceResult);

    return 0;
}