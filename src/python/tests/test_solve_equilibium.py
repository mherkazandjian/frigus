from __future__ import print_function
import numpy
from numpy.testing import assert_allclose
from astropy import units as u

from population import solveEquilibrium


def test_equilibrium_solver_snaity_2x2():

    M_matrix = numpy.array(
        [
            [0.2, 0.2],   # this row is not relevant since it is set to 1, 1
            [1.0, -1.0]   # in the solver
        ]
    )
    x = solveEquilibrium(M_matrix)

    expected_x_values = numpy.array([[0.5], [0.5]])
    assert_allclose(x, expected_x_values, rtol=1e-16, atol=0.0)



