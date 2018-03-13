# -*- coding: utf-8 -*-

#    linear.py is part of Frigus.

#    Frigus: software to compure the energy exchange in a multi-level system
#    Copyright (C) 2016-2018 Mher V. Kazandjian and Carla Maria Coppola

#    Frigus is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, version 3 of the License.    
#
#    Frigus is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Frigus.  If not, see <http://www.gnu.org/licenses/>.


"""module that implements helper functions for solving linear systems"""
import numpy
from numpy.linalg import solve, cond
import scipy

import mpmath
from mpmath import svd_r
mpmath.mp.dps = 50


def solve_linear_system_two_step(A, b, n_sub=1):
    """
    Solve the linaer system by sub-dividing into two linear systems

    The M-matrix property is exploited.

    :return: the solution of the linear system
    """
    n = n_sub

    C11 = A[0:n, 0:n]

    C12 = A[0:n, n:]

    A22 = A[n:, n:]
    A21 = A[n:, 0:n]

    C21 = scipy.linalg.solve(A22, A21)

    CD = C11 - numpy.dot(C12, C21)

    x1 = scipy.linalg.solve(CD, b[0:n])

    x2 = - numpy.dot(C21, x1.T)

    x = numpy.vstack((x1, x2))

    return x


def solve_equilibrium(m_matrix):
    """
    Solve for the equilibrium population densities given the right hand side of
    the linear system of the rate equations dx/dt as a matrix dx/dt = A.x
    where M_matrix = A in the case of this function. The first row of the A
    is replaces by the conservation equation.

    :param matrix_like m_matrix: The right hand side matrix of the rate
    equation as/home/carla an n x n matrix.
    :return: The population densities as a column vector.
    """

    sz = m_matrix.shape[0]

    #
    # set to zero rates that could be problematic while solving the system
    #
    # THRESHOLD_RATE = 1.0e-45
    # M_matrix[numpy.abs(M_matrix) < THRESHOLD_RATE] = 0.0

    # solving directly. replacing the first row with the conservation equation
    # i.e the sum of the independent variable is 1, i.e the sum of the
    # population levels is
    dxdt = numpy.zeros((sz, 1), 'f8')
    m_matrix[0, :], dxdt[0] = 1.0, 1.0

    # solving the system A.x = b
    # before solving, we will divide each row by the diagonal
    A, b = m_matrix, dxdt

    # ============ condition the linear system ========================
    #
    # scale the rows by normalizing w.r.t the diagonal element
    for i in numpy.arange(sz):
        A[i, :] = A[i, :] / A[i, i]

    #
    # DEBUG: examine the condition number of A
    #
    # cond = numpy.linalg.cond(A)
    # ============ done conditioning the linear system ===============

    try:
        x = solve_linear_system_two_step(A, b, n_sub=1)
    except scipy.linalg.LinAlgError as exc:
        print('solving the linear system with conditioning failed')
        print('due an singular matrix exception')
        print(exc)
        print('try to solve the system using extended precision')
        print('caution: this might take very long')
        x = solve_lu_mp(A, b)

    if x.any() < 0.0:
        print(
            'WARNING: found negative population densities\n'
            'accurary of the solution is not guaranteed.\n'
            'Check the linear system, rates, condition number..etc..\n'
        )

    return x


def solve_lu_mp(A, b):
    """
    Solve a linear system using mpmath

    :param ndarray A: The linear system
    :param ndarray b: The right hand side
    :return: ndarray
    """
    A_mp = mpmath.matrix([list(row) for row in A])
    b_mp = mpmath.matrix([list(row) for row in b])
    x_mp = mpmath.lu_solve(A_mp, b_mp)
    x = numpy.array(x_mp.tolist(), 'f8')

    return x


def solve_svd(A, b):
    """
    solve a linear system using svd decomposition using numpy

    :param ndarray A: The linear system
    :param ndarray b: The right hand side
    :return: ndarray
    """
    u, s, v = numpy.linalg.svd(A)

    c = numpy.dot(u.T, b)
    w = numpy.linalg.solve(numpy.diag(s), c)
    x = numpy.dot(v.T, w)

    return x


def solve_svd_mp(A, b):
    """
    solve a linear system using svd decomposition using mpmath

    :param ndarray A: The linear system
    :param ndarray b: The right hand side
    :return: ndarray
    """
    A_mp = mpmath.matrix([list(row) for row in A])
    b_mp = mpmath.matrix([list(row) for row in b])

    u, s, v = svd_r(A_mp)

    # x = V*((U'.b)./ diag(S))
    # x = V*(  c   ./ diag(S))
    # x = V*(       g        )

    c = u.T * b_mp
    w = mpmath.lu_solve(mpmath.diag(s), c)
    x_mp = v.T * w
    x = numpy.array(x_mp.tolist(), 'f8')

    return x
