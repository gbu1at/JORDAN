from functions import *
from main import *
import numpy as np
from smith_form import *
from sympy import Matrix

n = int(input())

TEST = []

for i in range(n):
    line = list(map(int, input().split()))
    TEST.append(line)

TEST = np.array(TEST)

n = TEST.__len__()

# CT = np.array(
#     [[randint(-10, 10) for _ in range(n)] for _ in range(n)]
# )
#
# CT_1 = np.linalg.inv(CT)

# TEST = np.array(
#     [[2, 0, 1],
#      [0, 6, 2],
#      [0, 0, 2]]
# )

# A = CT.dot(TEST).dot(CT_1)

A = TEST

n = A.__len__()

# A = TEST

smith_form = get_smith_form(A, n)

a = get_roots_in_smith_form(smith_form, n)

jung_diagrams = get_jung_diagrams(a)

# jung_diagrams = {-1 : {1 : 1}, 3 : {3 : 1}, 2 : {2 : 1}}

C = []

for root in jung_diagrams:
    jung_diagram = dict(sorted(jung_diagrams[root].items(), reverse=True))
    local_base_transition = get_jordan_cell(n, A, root, jung_diagram)
    C.extend(local_base_transition)

C = np.array(C).T
C_1 = np.linalg.inv(C)

print(C_1)
print(C)
print()
print()
result = C_1.dot(A).dot(C)

for i in range(n):
    print(*list(map(lambda x: round(x, 2), result[i])))
