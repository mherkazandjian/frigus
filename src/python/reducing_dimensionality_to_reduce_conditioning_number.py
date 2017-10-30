import numpy as np
import scipy.linalg as la

A = numpy.loadtxt('A.txt')
b = zeros((108, 1), 'f8')

C = np.copy(A)

# matrix well-conditioned removing the first 2 rows
# and columns; the original A matrix can then be written
# as:
# A = A11 A21
#     A21 A22
# The original matrix A can be multiplied by a matrix like:
# I 0
# 0 A22**{-1}
# (-A22 is a M-matrix, (-A22)**{-1}>=0)
# For this reason, the matrix is changed as -A
# Eventually, multiplying this matrix by A:
# I 0         *  A     = A11             A12
# 0 A22**{-1}            A22**{-1}*A21   I
# The linear system trasforms into:
# I 0         *  A * x =  I 0             1    = b*
# 0 A22**{-1}             0 A22**{-1}     0
#
# A11          A12  * x = b1
# A22**{-1}A21 I           0
#
#                            I
# A11             A12   * x1  = b1
# A22**{-1}*A21   I       x2    b2
#        C
# C = C11   C12
#     C21   I
# C11 x1 + C12 x2 = b1
# C21 x1 + x2 = 0
#
# x2 = -C21 x1
# C11 x1 -C12 C21 x1 = b1
# C11 x1 -C12 C21 x1 = b1
# (C11-C12 C21) x1 = b1
# ~1    ~1              x1  = 1
# 0     ~10^{-33}       x2    0

C = numpy.copy(A)

C21= la.solve(C[2:108, 2:108],C[2:108,0:1])

C11 = C[0:2,0:2]

C12 = C[0:2, 2:108]

CD = C11 - np.dot(C12, C21)

x1 = la.solve(CD, b[0:2])

x2 = - np.dot(C21, x1.T)

x1_list = x1[:][0].tolist()

x = [x1_copy[0][0], x1_copy[1][0]]

for i in arange(106):
    x2_only = x2_copy[i][0]
    x.append(x2_only)


