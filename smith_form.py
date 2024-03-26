import numpy as np
from sympy import factor, real_roots
from sympy.matrices.normalforms import smith_normal_form
from sympy import Poly, Matrix, QQ
from sympy.abc import x


def get_smith_form(A: np.array, n):
    M = Matrix([[A[i][j] if i != j else A[i][j] - x for j in range(n)] for i in range(n)])
    smith_matrix_form = smith_normal_form(M, QQ[x])
    return smith_matrix_form


def get_roots_in_smith_form(smith_form, n):
    a = []
    for i in range(n):
        v = []
        poly = smith_form.row(i).col(i)[0]
        if poly != 1:
            v = real_roots(poly)
        a.append(v)

    return a


def get_jung_diagrams(a):
    jung_diagrams_all_roots = {}

    for real_roots_in_smith in a:
        roots = set(real_roots_in_smith)
        for root in roots:
            cnt = real_roots_in_smith.count(root)
            if root in jung_diagrams_all_roots:
                if cnt in jung_diagrams_all_roots[root]:
                    jung_diagrams_all_roots[root][cnt] += 1
                else:
                    jung_diagrams_all_roots[root][cnt] = 1
            else:
                jung_diagrams_all_roots[root] = {cnt: 1}
    return jung_diagrams_all_roots

# v = {}
#
# n = 5
#
# A = Matrix([ [-4, 4, 0, 0, 0],
#              [0, 0, 0, 0, 0],
#              [0, 4, 4, 0, 0],
#              [3, -9, -4, 2, -1] ,
#              [1, 5, 4, 1, 4]])
#
# smith_form = get_smith_form(A, n)
# a = get_roots_in_smith_form(smith_form, n)
#
#
# print(get_jung_diagrams(a))
