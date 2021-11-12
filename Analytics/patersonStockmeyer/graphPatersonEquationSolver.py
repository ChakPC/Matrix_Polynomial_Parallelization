import os

os.system("rm output.txt")
for i in range(10, 1000, 100):
    print(f"Generating input matrix of size: {i}")
    os.system(f"python3 testGeneratorPatersonStockmeyer.py {i}")
    print(f"Running sylvester on matrix of size {i}")
    os.system("./paterson")
