#include "../complexMatrix.h"
#include "schurDecompositionSerial.h"

// Driver function for testing
int main()
{
    // Take Input
    int n;
    cin >> n;
    complexMatrix A(n, vector<complexNumber>(n, {0.0, 0.0}));
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            double a, b;
            cin >> a >> b;
            A[i][j] = {a, b};
        }
    }
    vector<complexMatrix> t = schurDecompositionSerial(A, 10);
    printf("\n Q: \n");
    printMatrix(t[0]);
    printf("\n Ak: \n");
    printMatrix(t[1]);
    return 0;
}