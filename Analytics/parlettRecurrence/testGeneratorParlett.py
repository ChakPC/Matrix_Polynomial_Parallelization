import random
import sys

sys.stdout = open("input.txt", "w")
SIZE_MATRIX = (int)(sys.argv[1])
COEFF_SIZE = 50

print(SIZE_MATRIX)
for i in range(SIZE_MATRIX):
    for j in range(SIZE_MATRIX):
        for _ in range(2):
            print(random.randint(1, 10), end=" ")
    print()

print(COEFF_SIZE)
for _ in range(COEFF_SIZE):
    for _ in range(3):
        print(random.randint(1, 10), end=" ")

sys.stdout.close()