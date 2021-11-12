import os

os.system("rm output.txt")
for i in range(10, 100, 10):
    print(f"Generating input matrix of size: {i}")
    os.system(f"python3 testGeneratorSylvesterEquationSolver.py {i}")
    print(f"Running sylvester on matrix of size {i}")
    os.system("./sylvester")
