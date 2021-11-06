// Header file for complex matrix operations

#include <bits/stdc++.h>
using namespace std;

#define complexNumber complex<double>
#define complexMatrix vector<vector<complexNumber>>

// + operaotr overloaded to add two complex matrices
inline complexMatrix operator+(complexMatrix &A, complexMatrix &B){
    int size = A.size();
    complexMatrix res(size, vector<complexNumber>(size, complexNumber(0.0, 0.0)));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            res[i][j] = A[i][j] + B[i][j];
        }
    }
    return res;
}

// * operator overloaded to multiply a complex number with a complex matrix
inline complexMatrix operator*(complexNumber coeff, complexMatrix &matrix){
    int size = matrix.size();
    complexMatrix res(size, vector<complexNumber>(size, {0.0, 0.0}));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            res[i][j] = coeff * matrix[i][j];
        }
    }
    return res;
}

// * operator overloaded to multiply two complex matrices
inline complexMatrix operator*(complexMatrix &A, complexMatrix &B){
    int size = A.size();
    complexMatrix res(size, vector<complexNumber>(size, {0.0, 0.0}));
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k < size; k++){
                res[i][j] = res[i][j] + A[i][k] * B[k][j];
            }
        }
    }
    return res;
}

// Function to print complex matrix
inline void printMatrix(complexMatrix &res){
    int n = res.size();
    int m = res[0].size();
    for(int i=0; i<n; i++){
        for(int j=0; j<m; j++){
            cout << "( " << real(res[i][j]) <<" ) + i( "<< imag(res[i][j]) <<" )  ";
        }
        cout << "\n";
    }
}

// Function to print complex number
inline void printComplex(complexNumber res) {
    cout << "(" << real(res) << " + i" << imag(res) << ") ";
}

// Function to return transpose of a complex matrix
// Returns the conjugate transpose
inline complexMatrix transposeMatrix(complexMatrix& inputMatrix) {
    int n = inputMatrix.size();
    complexMatrix newMatrix(n, vector<complexNumber>(n, complexNumber(0, 0)));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double Re = real(inputMatrix[j][i]);
            double Im = imag(inputMatrix[j][i]);
            newMatrix[i][j] = complexNumber(Re, -1.0 * Im);
        }
    }
    return newMatrix;
}

inline complexMatrix identityMatrix(int n) {
    complexMatrix res(n, vector<complexNumber>(n, complexNumber(0.0, 0.0)));
    for(int i=0; i<n; i++){
        res[i][i] = complexNumber(1.0, 0.0);
    }
    return res;
}

// Function to compute kronecker product of A and B
inline complexMatrix computeKroneckerProduct(complexMatrix &A, complexMatrix &B){
    int rA = A.size(), cA = A[0].size();
    int rB = B.size(), cB = B[0].size();
    complexMatrix res(rA * rB, vector<complexNumber>(cA * cB, complexNumber(0, 0)));
    for(int i=0; i<rA; i++){
        for(int j=0; j<cA; j++){
            complexMatrix Aij_B = A[i][j] * B;
            int rowStart = i * rB, rowEnd = (i + 1) * rB;
            int colStart = j * cB, colEnd = (j + 1) * cB;
            for(int x = rowStart; x < rowEnd; x++){
                for(int y = colStart; y < colEnd; y++){
                    res[x][y] = Aij_B[x - rowStart][y - colStart];
                }
            }
        }
    }
    return res;
}

// Function to stack columns of a matrix atop each other from left to right
inline complexMatrix columnStack(complexMatrix &inputMatrix){
    int p = inputMatrix.size(), q = inputMatrix[0].size();
    complexMatrix res(p * q, vector<complexNumber>(1, complexNumber(0, 0)));
    int pos = p * q - 1;
    for(int i = 0; i < q; i++){
        for(int j = p-1; j >= 0; j--){
            res[pos][0] = inputMatrix[j][i];
            pos--;
        }
    }
    return res;
}

// Function to convert column stacked matrix back to standard form
inline complexMatrix reverseStacking(complexMatrix &inputMatrix, int rows){
    int prod = inputMatrix.size();
    int cols = prod / rows;
    complexMatrix res(rows, vector<complexNumber>(cols, complexNumber(0, 0)));
    int pos = 0;
    for(int i = cols-1; i>=0; i--){
        for(int j = 0; j < rows; j++){
            res[j][i] = inputMatrix[pos][0];
            pos++;
        }
    }
    return res;
}

// Function to get cofactor of A[p][q] in temp[][]. n is current dimension of A[][]
inline complexMatrix getCofactor(complexMatrix &A, int p, int q, int n){
    complexMatrix res(A.size(), vector<complexNumber>(A.size(), complexNumber(0, 0)));
	int i = 0, j = 0;
	for (int row = 0; row < n; row++){
		for (int col = 0; col < n; col++){
			if (row != p && col != q){
				res[i][j] = A[row][col];
                j++;
				if (j == n - 1){
					j = 0;
					i++;
				}
			}
		}
	}
    return res;
}

// Recursive function for finding determinant of matrix.
// n is current dimension of A
inline complexNumber determinant(complexMatrix A, int n){
	if (n == 1) return A[0][0];
	complexNumber D(0, 0);
    complexMatrix temp(A.size(), vector<complexNumber>(A.size(), complexNumber(0, 0)));
	int sign = 1;
	for (int f = 0; f < n; f++){
		temp = getCofactor(A, 0, f, n);
        complexNumber num = complexNumber(sign, 0) * A[0][f];
		D += num * determinant(temp, n - 1);
		sign = -sign;
	}
	return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
inline complexMatrix adjoint(complexMatrix &A){
    complexMatrix adj(A.size(), vector<complexNumber>(A.size(), complexNumber(0, 0)));
	if (A.size() == 1){
		adj[0][0] = complexNumber(1, 0);
		return adj;
	}
    complexMatrix temp(A.size(), vector<complexNumber>(A.size(), complexNumber(0, 0)));
	int sign = 1;
	for (int i=0; i<A.size(); i++){
		for (int j=0; j<A.size(); j++){
			temp = getCofactor(A, i, j, A.size());
			sign = ((i+j)%2==0)? 1: -1;
			adj[j][i] = complexNumber(sign, 0)*(determinant(temp, A.size()-1));
		}
	}
    return adj;
}

// Function to calculate and store inverse, returns false if
// matrix is singular
inline complexMatrix inverse(complexMatrix A)
{
    complexMatrix inverse(A.size(), vector<complexNumber>(A.size(), complexNumber(0, 0)));
	complexNumber det = determinant(A, A.size());
    cout << "\n Det: " << real(det) <<" " << imag(det) <<"\n";
	complexMatrix adj = adjoint(A);
	for (int i=0; i<A.size(); i++)
		for (int j=0; j<A.size(); j++)
			inverse[i][j] = adj[i][j] / det;
	return adj;
}