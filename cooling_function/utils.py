def linear_2d_index(i, j):
    """
    given two integer array returns an array of the same size mapping the
    elements to a unique integer (assuming there are non repetitive i,j pairs.

    :param i: an integer array
    :param j: an integer array
    :return: array give i,j pair unique integer label
    """
    return j*(i.max() + 1) + i
