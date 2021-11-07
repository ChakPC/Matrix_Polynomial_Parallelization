#include "complexMatrix.h"
#include "Paterson_Stockmeyer/patersonStockmeyerSerial.h"
#include "Schur_Decomposition/schurDecompositionSerial.h"
#include "Sylvester_Equation_Solver/sylvesterEquationSolverSerial.h"

// Function to compute matrix polynomial using Parlett recurrence
complexMatrix parlettRecurrenceSerial(complexMatrix &A, vector<complexNumber> &coeff, int numOfIterations, int PS_p, int PS_s){
    // Schur Decomposition
    vector<complexMatrix> schur = schurDecompositionSerial(A, numOfIterations);
    complexMatrix Q = schur[0];
    complexMatrix T = schur[1];

    // Compute function for diagonal elements
    int size = Q.size();
    for(int i=0; i<size; i++){
        if(i + 1 < size && abs(real(Q[i+1][i])) > 1e-15 && abs(real(Q[i+1][i])) > 1e-15){

        }
    }
}

int main(){
    // Take Input
    int n;
    cin >> n;
    complexMatrix A(n, vector<complexNumber>(n, {0.0, 0.0}));
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    int d;
    cin >> d;
    vector<complexNumber> coeff(d+1, {0.0, 0.0});
    for(int i=0; i<=d; i++){
        double a, b;
        cin >> a >> b;
        coeff[i] = {a, b};
    }

    // Output from Paterson Stockmeyer
    complexMatrix patersonStockmeyerResult = patersonStockmeyerSerial(A, coeff, sqrt(d)+1, sqrt(d)+1);
    printMatrix(patersonStockmeyerResult);

    // Output from Parlett Recurrence
    complexMatrix parlettRecurrenceResult = parlettRecurrenceSerial(A, coeff, 10, sqrt(d)+1, sqrt(d)+1);
    printMatrix(parlettRecurrenceResult);

    return 0;
}
