#include "complexMatrix.h"
#include "Paterson_Stockmeyer/patersonStockmeyerSerial.h"
#include "Schur_Decomposition/schurDecompositionSerial.h"
#include "Sylvester_Equation_Solver/sylvesterEquationSolverSerial.h"

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
    

    return 0;
}
