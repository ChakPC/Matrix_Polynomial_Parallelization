<h1 align="center">Matrix Polynomial Parallelization</h1>

Optimized and Parallelized Code for Solving Matrix Polynomials  


### Problem Statement
Solving the given equation
<p align="center">
  <img src="https://user-images.githubusercontent.com/59741135/142240642-d246624d-8940-464a-b9a3-0d3f548a04e4.png" style="filter: invert(1);"/>
</p>


## Installation


```bash
git clone https://github.com/ChakPC/Matrix_Polynomial_Parallelization.git
cd Matrix_Polynomial_Parallelization
```

## Testing


```bash
cd Testing
g++ patersonStockmeyer.cpp -fopenmp             # compiling paterson stockmeyer with openmp
./a.out                                         # executing compiled patersonStockmeyer
g++ schurDecomposition.cpp -fopenmp             # compiling schur decomposition with openmp
./a.out                                         # executing schur decompositions (may take some time to execute)
g++ sylvesterEquationSolver.cpp -fopenmp        # compiling sylvester with openmp
./a.out                                         # executing sylvester (may take a few minutes for serial calculations)
g++ parletRecurrence.cpp -fopenmp               # compiling parlett Recurrence with fopenmp
./a.out                                         # executing parlett (may take a few minutes for serial calculations)
```

