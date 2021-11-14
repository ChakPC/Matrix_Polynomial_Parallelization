<h1 align="center">Matrix Polynomial Parallelization</h1>

Optimized and Parallelized Code for Solving Matrix Polynomials  


### Problem Statement
Solving the given equation
<p align="center">
  <img src="https://latex.codecogs.com/svg.latex?q%28x%29%20%3D%20%5Csum%5Climits_%7Bi%3D1%7D%5Ed%20c_ix%5Ed%20%3D%20c_0%20+%20c_1x%20+%20c_2x%5E2%20+%20...c_dx%5Ed" style="filter: invert(1);"/>
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

