from __future__ import print_function
import numpy
from numpy.testing import assert_allclose

from frigus.population import solve_equilibrium


def test_equilibrium_solver_snaity_2x2():

    m_matrix = numpy.array(
        [
            [0.2, 0.2],   # this row is not relevant since it is set to 1, 1
            [1.0, -1.0]   # in the solver
        ]
    )
    x = solve_equilibrium(m_matrix)

    expected_x_values = numpy.array([[0.5], [0.5]])
    assert_allclose(x, expected_x_values, rtol=1e-15, atol=0.0)


def test_equilibrium_solver_snaity_3x3():

    m_matrix = numpy.array(
        [
            [0.2, 0.2, 0.2],   # this row is not relevant since it is set to
            [1.0, -1.0, 1.0],  # 1, 1, 1 in the solver
            [3.0, -2.0, -9.0]
        ]
    )
    x = solve_equilibrium(m_matrix)

    expected_x_values = numpy.array([[11.0/24.0], [1.0/2.0], [1.0/24.0]])
    assert_allclose(x, expected_x_values, rtol=1e-14, atol=0.0)



