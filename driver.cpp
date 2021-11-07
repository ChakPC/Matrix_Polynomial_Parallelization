#include "complexMatrix.h"
#include "Paterson_Stockmeyer/patersonStockmeyerSerial.h"
#include "Schur_Decomposition/schurDecompositionSerial.h"
#include "Sylvester_Equation_Solver/sylvesterEquationSolverSerial.h"


pair<int,map<pair<int,int> , complexMatrix>> constructBlocks(complexMatrix &inputMatrix){
    map <pair<int,int>,complexMatrix> res;
    map<pair<int,int>, pair<int,int>> limits;
    int size = inputMatrix.size();
    int count=0;
    for(int i=0; i<size; i++){
        count++;
        if(i + 1 < size && (abs(real(inputMatrix[i+1][i])) > 1e-15 || abs(imag(inputMatrix[i+1][i])) > 1e-15)){
            complexMatrix block(2,vector<complexNumber>(2));
            block[0][0] = inputMatrix[i][i];
            block[0][1] = inputMatrix[i][i+1];
            block[1][0] = inputMatrix[i+1][i];
            block[1][1] = inputMatrix[i+1][i+1];
            res[{count,count}] = block;
            limits[{count,count}] = {i,i};
            i++;
        }
        else{
            complexMatrix block(1,vector<complexNumber>(1));
            block[0][0] = inputMatrix[i][i];
            res[{count,count}] = block;
            limits[{count,count}] = {i,i};
        }
    }

    for(int i=1;i<=count;i++){
        for(int j=i+1;j<=count;j++){
            int row = res[{i,i}].size();
            int col = res[{j,j}].size();
            complexMatrix block(row,vector<complexNumber>(col));
            int rowStart = limits[{i,i}].first;
            int colStart = limits[{j,j}].second;
            for(int x=rowStart; x<rowStart+row; x++){
                for(int y=colStart; y<colStart+col; y++){
                    block[x-rowStart][y-colStart] = inputMatrix[x][y];
                }
            }
            res[{i,j}] = block;
        }
    }
    return {count,res};
}

// Function to compute matrix polynomial using Parlett recurrence
complexMatrix parlettRecurrenceSerial(complexMatrix &A, vector<complexNumber> &coeff, int numOfIterations, int PS_p, int PS_s){
    // Schur Decomposition
    vector<complexMatrix> schur = schurDecompositionSerial(A, numOfIterations);
    complexMatrix Q = schur[0];
    complexMatrix T = schur[1];

    // Compute function for diagonal elements
    pair<int,map<pair<int,int> , complexMatrix>> res = constructBlocks(T);
    int numOfBlocks = res.first;
    map<pair<int,int> , complexMatrix> blocks = res.second;
    int size = A.size();
    complexMatrix functionOnT(size,vector<complexNumber>(size,complexNumber(0,0)));
    for(int i=1;i<=numOfBlocks;i++){
        complexMatrix computedBlock = patersonStockmeyerSerial(blocks[{i,i}], coeff, PS_p, PS_s);
        if(computedBlock.size()){
            
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

    // for(int i=0;i<4;i++){
    //     for(int j=i;j<4;j++){
    //         cout << "\n( " << i << " " << j << " )" << endl;
    //         printMatrix(temp[{i,j}]);
    //     }
    // }

    // int d;
    // cin >> d;
    // vector<complexNumber> coeff(d+1, {0.0, 0.0});
    // for(int i=0; i<=d; i++){
    //     double a, b;
    //     cin >> a >> b;
    //     coeff[i] = {a, b};
    // }



    // Output from Paterson Stockmeyer
    // complexMatrix patersonStockmeyerResult = patersonStockmeyerSerial(A, coeff, sqrt(d)+1, sqrt(d)+1);
    // printMatrix(patersonStockmeyerResult);

    // // Output from Parlett Recurrence
    // complexMatrix parlettRecurrenceResult = parlettRecurrenceSerial(A, coeff, 10, sqrt(d)+1, sqrt(d)+1);
    // printMatrix(parlettRecurrenceResult);

    return 0;
}
