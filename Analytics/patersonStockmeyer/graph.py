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
for i in range(1, 1000, 100):
    X.append(i)

plt.figure(figsize=(9 * 1.5, 6 * 1.5))
plt.plot(X, P, color="blue", label="serial")
plt.plot(X, S, color="green", label="parallel")
plt.legend()
plt.xlabel("Matrix Dimension")
plt.ylabel("Time in seconds")
plt.title("Paterson Stockmeyer")
plt.show()