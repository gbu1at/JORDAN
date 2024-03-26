import numpy
import numpy as np
from random import randint

Matrix = np.array
Vector = np.array
SPREAD = 100
eps = 1e-6


def pow(B: Matrix, n):
    if n == 0:
        sz = B.__len__()
        return np.eye(sz)

    if n % 2 == 0:
        A = pow(B, n // 2)
        return A.dot(A)

    return B.dot(pow(B, n - 1))


def get_random_vector(n: int, base_random_vector: np.array):
    """
        n - размерность вектора
    """

    vector = np.array([0] * n)

    for i in range(base_random_vector.__len__()):
        vector = vector + randint(-100, 100) * base_random_vector[i]

    return vector


def is_linear_dependence_vectors(A: Matrix) -> bool:
    """
        возвращает true если вектора линейно зависимы
        иначе false
    """

    a = A.copy()

    n = len(a)
    m = len(A[0])

    for col in range(m):
        i = -1
        for row in range(col, n):
            if abs(a[row][col]) > eps:
                i = row

        for row in range(col, n):
            if row != i:
                k = a[row][col]
                a[row] *= a[i][col]
                a[row] -= a[i] * k

        if i != -1:
            a[col], a[i] = a[i].copy(), a[col].copy()
        else:
            return True

    return False


def set_linear_dependence_coefficients(A: Matrix) -> list:
    """
        вернет набор коэффициентов в линейной зависимости векторов
    """
    a = A.copy()

    n = len(a)
    m = len(A[0])

    for col in range(m):
        i = -1
        for row in range(col, n):
            if abs(a[row][col]) > eps:
                i = row

        if i != -1:
            for row in range(col, n):
                if row != i:
                    k = a[row][col]
                    a[row] *= a[i][col]
                    a[row] -= a[i] * k

            a[col], a[i] = a[i].copy(), a[col].copy()

    coeff = [0] * m

    for i in range(m - 1, -1, -1):
        s = 0
        for j in range(i + 1, m):
            s += a[i][j] * coeff[j]
        if abs(a[i][i]) < eps:
            coeff[i] = 1
        else:
            coeff[i] = -s / a[i][i]

    return coeff


def is_null_vector(v: Vector) -> bool:
    """
        вернет true если вектор нулевой
        иначе false
    """
    n = len(v)

    for i in range(n):
        if abs(v[i]) > eps:
            return False
    return True


def get_kernel_linear_operator(B):
    a = B.copy()

    n = len(a)

    for col in range(n):
        i = -1
        for row in range(col, n):
            if abs(a[row][col]) > eps:
                i = row


        if i != -1:
            for row in range(col, n):
                if row != i:
                    k = a[row][col]
                    a[row] *= a[i][col]
                    a[row] -= a[i] * k
            a[col], a[i] = a[i].copy(), a[col].copy()

    idxes = []

    for i in range(n):
        if abs(a[i][i]) < eps:
            idxes.append(i)

    rank = idxes.__len__()

    coeff = [[0 for _ in range(n)] for _ in range(rank)]

    for _ in range(rank):
        coeff[_][idxes[_]] = 1
        for i in range(n - 1, -1, -1):
            s = 0
            for j in range(i + 1, n):
                s += a[i][j] * coeff[_][j]
            if idxes[_] != i:
                if abs(a[i][i]) < eps:
                    coeff[_][i] = 0
                else:
                    coeff[_][i] = -s / a[i][i]

    return np.array(coeff)


def get_base_for_operator(B: Matrix, k):
    B_ = pow(B, k)
    return get_kernel_linear_operator(B_)
