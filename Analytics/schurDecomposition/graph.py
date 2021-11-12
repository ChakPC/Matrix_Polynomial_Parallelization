import matplotlib.pyplot as plt
from numpy import double

with open('output.txt') as f:
    lines = f.readlines()

S = []
P = []

for i in range(len(lines)):
    a = lines[i].split(":")[1]
    if (i % 2):
        S.append(double(a))
    else:
        P.append(double(a))

X = []
for i in range(10, 80, 10):
    X.append(i)

plt.figure(figsize=(9 * 1.5, 6 * 1.5))
plt.plot(X, P, color="blue", label="parallel")
plt.plot(X, S, color="green", label="serial")
plt.xlabel("Matrix Dimension")
plt.ylabel("Time in seconds")
plt.title("Schur Decomposition")
plt.show()