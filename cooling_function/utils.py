import numpy
import os
from numpy.testing import assert_approx_equal

def linear_2d_index(i, j, n_i=None):
    """
    given two integer array returns an array of the same size mapping the
    elements to a unique integer (assuming there are non repetitive i,j pairs.
    it is assumed that the array i and j have a range between 0 and i.max() + 1
    and j.max() + 1 respectively.

    :param i: an integer array
    :param j: an integer array
    :param int n_i: manually set the "number" of i element along the i axis
     otherwise the i.max() + 1 is used. (0 indexing is used)
    :return: array give i,j pair unique integer label
    """

    if n_i is None:
        n_i = i.max() + 1

    return j*n_i + i


def find_matching_indices(v1, v2, check=True):
    """find indices in v1 of intersecting elements of v2 in v1. It is assumed
     that all the entries is v2 are a subset of v1. The elements of v2 must
     be unique. The returned indices array
     (r) is the same size as v2, is such that :

     v1[r[i]] = v2[i]

     in other words :

     v1[r] == v2

     or

     v2 == v1[r]

    :param array v1: reference array
    :param array v2: array that will be looked up in v1
    :param bool check: checks whether the unique elements of v2 are a subset
     of v1.

    .. see-also:: tests/test_find_matching_indices.py

    .. see-also:: numpy.in1d, numpy.intersect1d, numpy.setdiff1d

    .. code-block:: python

        ind1   = [0, 1, 2, 3, 4, 5,  6]
        v1     = [5, 1, 3, 7, 6, 8, -2]
        ind1s  = [6, 1, 2, 0, 3, 4, 5]  # indices in v1

        v2     = [7, -2, 3, 5, 5, 5, 1, 1, 7, -2, 6, 7, 8]
        v2us   = [-2, 1, 3, 5, 7, 6, 8]  # sorted unique v2

        # returned array
        print find_matching_indices(v1, v2)
        >>>  [3,  6, 2, 0, 0, 0, 1, 1, 3,  6, 4, 3, 5]

        # now if we want to change the replace the 7 with a -9 in v1 and do the
        # the same for v2, we can achieved this via:
        v1[ v1 == 7 ] = -9
        v2[:] = v1[r]

        print v1
        print v2

        # note that here all the values of v1[r] are assigned to v2, i.e all the
        # changes other than replacing 7 with -9 would be transferred to v2

        # ------------------------------------------------------------------
        # a second example where not all the values of v1 are used in v2
        v1   = [5, 1, 3, 7, 6, 8, -2]
        ind1 = [0, 1, 2, 3, 4, 5,  6]

        v2 = [7, -2, 3, 1, 1, 7, -2, 6, 7, 8]

        v2us  = [-2, 1, 3, 7, 6, 8]  # unique sorted v2
        ind1s = [ 6, 1, 2, 3, 4, 5]   # indices in v1

        # returned array
        print find_matching_indices(v1, v2)
        >>>  [3, 6, 2, 1, 1, 3, 6, 4, 3, 5]

        # now if we want to change the replace the 7 with a -9 in v1 and do the
        # the same for v2, we can achieved this via:
        v1[ v1 == 7 ] = -9
        v2[:] = v1[r]

        print v1
        print v2

        # note that here all the values of v1[r] are assigned to v2, i.e all the
        # changes other than replacing 7 with -9 would be transferred to v2
    """

    # get the unique values in v1 along with the indices array that constructs
    # v1u from v1, i.e v1u = v1[inds_rec_v1]
    v1u, inds_rec_v1 = numpy.unique(v1, return_index=True)
    assert v1u.size == v1.size

    # getting the unique elements of v2 and the array of indices to recover it
    # from the unique values
    v2u, inds_rec_v2_from_v2u = numpy.unique(v2, return_inverse=True)

    # indices of v1u that have matches in v1
    inds_v1u_of_v2u_in_v1u = numpy.where(numpy.in1d(v1u, v2u))[0]

    # indices in v1 unique (that is also sorted) of elements in v2 with matches
    # in v1
    inds_v1_of_v2u_in_v1 = inds_v1u_of_v2u_in_v1u[inds_rec_v2_from_v2u]
    # assert v1u[inds_v1_of_v2u_in_v1] == v2   # DEBUG

    if check is True:
        if inds_v1_of_v2u_in_v1.size != v2.size:
            raise ValueError("not all unique entries in v2 are in v1")

    # indices in v1 of elements in v2 with matches in v1 v1[r] = v2
    inds_v1_of_v2_in_v1 = inds_rec_v1[inds_v1_of_v2u_in_v1]

    return inds_v1_of_v2_in_v1


def load_ascii_matrix_data():
    A_matrix_mher = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/mher/A_matrix.txt'))
    B_matrix_mher = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/mher/B_matrix.txt'))
    J_nu_matrix_mher = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/mher/J_nu_matrix.txt'))
    K_matrix_mher = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/mher/K_nc_matrix.txt'))
    O_matrix_mher = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/mher/O_matrix.txt'))
    D_matrix_mher = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/mher/D_matrix.txt'))
    M_matrix_mher = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/mher/M_matrix.txt'))


    A_matrix_carla = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/carla/A_matrix.txt')).T
    B_matrix_carla = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/carla/B_matrix.txt')).T
    J_nu_matrix_carla = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/carla/J_nu_matrix.txt')).T
    K_matrix_carla = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/carla/K_nc_matrix.txt')).T
    O_matrix_carla = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/carla/O_matrix.txt')).T
    D_matrix_carla = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/carla/D_matrix.txt')).T
    M_matrix_carla = numpy.loadtxt(
        os.path.expanduser(
            '~/dropbox/Dropbox/us/cooling_function/carla/M_matrix.txt')).T


    rel_diff_A = 1.0 - A_matrix_carla / A_matrix_mher
    rel_diff_A[numpy.where(numpy.isnan(rel_diff_A))] = 0.0
    assert rel_diff_A.sum() < 1e-14

    rel_diff_B = 1.0 - B_matrix_carla / B_matrix_mher
    rel_diff_B[numpy.where(numpy.isnan(rel_diff_B))] = 0.0
    assert numpy.fabs(rel_diff_B.sum()) < 1e-4
    assert numpy.fabs(rel_diff_B.max()) < 1e-7

    rel_diff_J_nu = 1.0 - J_nu_matrix_carla / J_nu_matrix_mher
    rel_diff_J_nu[numpy.where(numpy.isnan(rel_diff_J_nu))] = 0.0
    assert numpy.fabs(rel_diff_J_nu).max() < 1e-6

    rel_diff_K = 1.0 - K_matrix_carla / K_matrix_mher
    rel_diff_K[numpy.where(numpy.isnan(rel_diff_K))] = 0.0
    assert numpy.fabs(rel_diff_K).max() < 1e-12

    rel_diff_O = 1.0 - O_matrix_carla / O_matrix_mher
    rel_diff_O[numpy.where(numpy.isnan(rel_diff_O))] = 0.0
    assert numpy.fabs(rel_diff_O).max() < 1e-12

    rel_diff_D = 1.0 - D_matrix_carla / D_matrix_mher
    rel_diff_D[numpy.where(numpy.isnan(rel_diff_D))] = 0.0
    assert numpy.fabs(rel_diff_D).max() < 1e-12

    rel_diff_M = 1.0 - M_matrix_carla / M_matrix_mher
    rel_diff_M[numpy.where(numpy.isnan(rel_diff_M))] = 0.0
    assert numpy.fabs(rel_diff_M).max() < 1e-12

    print('done checking')