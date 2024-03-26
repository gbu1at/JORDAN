import numpy as np



M = [[-4, 4, 0, 0, 0],
     [0, 0, 0, 0, 0],
     [0, 4, 4, 0, 0],
     [3, -9, -4, 2, -1] ,
     [1, 5, 4, 1, 4]]


M = np.array(M)



# print(M.dot(M))
#
#
# A = [[0, 4, 0, 0, 0],
#      [0, 4, 0, 0, 0],
#      [0, 4, 8, 0, 0],
#      [3, -9, -4, 6, -1] ,
#      [1, 5, 4, 1, 8]]
#
# A = np.array(A)
#
# v = [[0],
#      [0],
#      [0],
#      [1]]
#
# v = np.array(v)
#
#
#
# C = [[], [], [], []]
# C[0] = v
# C[1] = A.dot(C[0])
# C[2] = A.dot(C[1])
# C[3] = np.array([[-1], [-5], [0], [2]])

A = [[-7, 4, 0, 0, 0],
     [0, -3, 0, 0, 0],
     [0, 4, 1, 0, 0],
     [3, -9, -4, -1, -1] ,
     [1, 5, 4, 1, 1]]

A = np.array(A)


print(A.dot(np.array([0, 0, 0, -7, 3])))

C = [[-49, 0, 0, 25, 3],
     [1, 1, -1, 2, -2],
     [0, 0, 1, -4, 4],
     [0, 0, 0, -7, 3],
     [0, 0, 0, 4, -4]]


C = np.array(C).T
# print(C)

C_1 = np.linalg.inv(C)


print(C_1)

print(C_1.dot(M).dot(C))