from itertools import chain
from functions import *
import numpy as np

Matrix = np.array


def operator(A: Matrix, l) -> Matrix:
    n = len(A)
    B = A.copy()
    for i in range(n):
        B[i][i] = A[i][i] - l

    return B


def get_jordan_cell(n: int, A: Matrix, l, jung_diagram: dict):
    """
            A - изеачальная мтарица
            l - осбственное число
            jung_diagram - диаграмма юнга для этого собственного числа
            эта диаграмма представляет собой словаь из длины соотвествуюшего блока и кол-ва этого блока в клетке
    """

    B = operator(A, l)

    # необзодимо найти базис в котором матрица B будет иметь якобы лильпотентный
    # вид то есть она будет занулять все вектора из простанства порожденого этоим базисом

    # В принципе алгоритм отличается от нильпотентных матриц только выбором рандомного базиса
    # в остальном все остается точно такое же

    all_base = [[]]

    # мы создали первый набор векторов максимального порядка
    # далее нам необходимо проверуть ту же стратегию для остальных векторов

    it_last_base = 0

    max_deg = list(jung_diagram.items())[0][0]

    base_random_vector = get_base_for_operator(B, max_deg)

    for rank_local_base, count_local_base in jung_diagram.items():

        for _ in range(count_local_base):

            u = get_random_vector(n, base_random_vector)

            v = u.copy()

            for _ in range(rank_local_base):
                v = B.dot(v)

            # v - это вектор из которого мы зотим вытащить коэффициенты базисных векторов. чтобы он обнулялся

            T = [v]

            for i in range(all_base.__len__()):
                for j in range(rank_local_base, all_base[i].__len__()):
                    T.append(all_base[i][j])
            #                     добавили векто. который может присутстовать в разложении векторва Bv

            T = np.array(T).transpose()
            assert is_linear_dependence_vectors(T)

            coeff = set_linear_dependence_coefficients(T)

            it_coeff = 0
            u = coeff[it_coeff] * u

            q = coeff[it_coeff] * v

            for i in range(all_base.__len__()):
                for j in range(rank_local_base, all_base[i].__len__()):
                    it_coeff += 1
                    u = u + coeff[it_coeff] * all_base[i][j - rank_local_base]
                    q = q + coeff[it_coeff] * all_base[i][j]

            all_base.append([])
            for _ in range(rank_local_base):
                all_base[it_last_base].append(u)
                u = B.dot(u)

            it_last_base += 1

            assert is_null_vector(u)

    all_base_ = []
    for line in all_base:
        for vec in line:
            all_base_.append(vec)
    all_base = all_base_
    all_base = np.array(all_base).T

    return all_base.T
