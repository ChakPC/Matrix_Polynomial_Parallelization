import random
import sys

sys.stdout = open("input.txt", "w")
SIZE_MATRIX = 30

print(SIZE_MATRIX)
for i in range(SIZE_MATRIX):
    for j in range(SIZE_MATRIX):
        for _ in range(2):
            print(random.randint(1, 10), end=" ")
    print()

print(SIZE_MATRIX)
for i in range(SIZE_MATRIX):
    for j in range(SIZE_MATRIX):
        for _ in range(2):
            print(random.randint(1, 10), end=" ")
    print()

# print(SIZE_MATRIX)
for i in range(SIZE_MATRIX):
    for j in range(SIZE_MATRIX):
        for _ in range(2):
            print(random.randint(1, 10), end=" ")
    print()

sys.stdout.close()