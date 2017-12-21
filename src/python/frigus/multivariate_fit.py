import numpy as np
from scipy.optimize import curve_fit
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def func(X, D00,  D10, D20, D30, D40, \
            D01,  D11, D21, D31, D41, \
            D02,  D12, D22, D32, D42, \
            D03,  D13, D23, D33, D43, \
            D04,  D14, D24, D34, D44):
    T, nh = X

    ####
    lin1 = D00
    lin2 = D10 * T
    lin3 = D20 * T**2
    lin4 = D30 * T**3
    lin5 = D40 * T**4
    ####
    ####
    lin6 = D01 * nh
    lin7 = D02 * nh**2
    lin8 = D03 * nh**3
    lin9 = D04 * nh**4
    ####
    ####
    lin10 = D11 * T * nh
    lin11 = D21 * (T**2) * nh
    lin12 = D31 * (T**3) * nh
    lin13 = D41 * (T**4) * nh
    ####
    ####
    lin14 = D12 * T * (nh**2)
    lin15 = D22 * (T**2) * (nh**2)
    lin16 = D32 * (T**3) * (nh**2)
    lin17 = D42 * (T**4) * (nh**2)
    ####
    ####
    lin18 = D13 * T * (nh**3)
    lin19 = D23 * (T**2) * (nh**3)
    lin20 = D33 * (T**3) * (nh**3)
    lin21 = D43 * (T**4) * (nh**3)
    ####
    ####
    lin22 = D14 * T * (nh**4)
    lin23 = D24 * (T**2) * (nh**4)
    lin24 = D34 * (T**3) * (nh**4)
    lin25 = D44 * (T**4) * (nh**4)
    ####

    sum = (lin1 + lin2 + lin3 + lin4 + lin5 +
          lin6 + lin7 + lin8 + lin9 + lin10 +
          lin11 + lin12 + lin13 + lin14 + lin15 +
          lin16 + lin17 + lin18 + lin19 + lin20 +
          lin21 + lin22 + lin23 + lin24 + lin25)
    return sum


def fit_lambda(x_fit, y_fit, data_to_fit):
    # D00 = - 42.57688
    # D10 = + 21.93385 #* lt_kin ** 1 * ln_hd ** 0
    # D20 = - 10.19097 #* lt_kin ** 2 * ln_hd ** 0
    # D30 = + 2.19906 #* lt_kin ** 3 * ln_hd ** 0
    # D40 = - 0.17334 #* lt_kin ** 4 * ln_hd ** 0
    # D01 = + 0.92433 #* lt_kin ** 0 * ln_hd ** 1
    # D11 = + 0.77952 #* lt_kin ** 1 * ln_hd ** 1
    # D21 = - 0.54263 #* lt_kin ** 2 * ln_hd ** 1
    # D31 = + 0.11711 #* lt_kin ** 3 * ln_hd ** 1
    # D41 = - 0.00835 #* lt_kin ** 4 * ln_hd ** 1
    # D02 = + 0.54962 #* lt_kin ** 0 * ln_hd ** 2
    # D12 = - 1.06447 #* lt_kin ** 1 * ln_hd ** 2
    # D22 = + 0.62343 #* lt_kin ** 2 * ln_hd ** 2
    # D32 = - 0.13768 #* lt_kin ** 3 * ln_hd ** 2
    # D42 = + 0.0106 #* lt_kin ** 4 * ln_hd ** 2
    # D03 = - 0.07676 #* lt_kin ** 0 * ln_hd ** 3
    # D13 = + 0.11864 #* lt_kin ** 1 * ln_hd ** 3
    # D23 = - 0.07366 #* lt_kin ** 2 * ln_hd ** 3
    # D33 = + 0.01759 #* lt_kin ** 3 * ln_hd ** 3
    # D43 = - 0.001482 #* lt_kin ** 4 * ln_hd ** 3
    # D04 = + 0.00275 #* lt_kin ** 0 * ln_hd ** 4
    # D14 = - 0.00366 #* lt_kin ** 1 * ln_hd ** 4
    # D24 = + 0.002514 #* lt_kin ** 2 * ln_hd ** 4
    # D34 = - 0.000666317 #* lt_kin ** 3 * ln_hd ** 4
    # D44 = + 0.000061926 #* lt_kin ** 4 * ln_hd ** 4)

    D00 = - 42.
    D10 = + 10. #* lt_kin ** 1 * ln_hd ** 0
    D20 = - 10. #* lt_kin ** 2 * ln_hd ** 0
    D30 = + 1. #* lt_kin ** 3 * ln_hd ** 0
    D40 = - 0.1 #* lt_kin ** 4 * ln_hd ** 0

    D01 = + 1.
    D11 = + 1. #* lt_kin ** 1 * ln_hd ** 1
    D21 = - 0.1 #* lt_kin ** 2 * ln_hd ** 1
    D31 = + 0.1 #* lt_kin ** 3 * ln_hd ** 1
    D41 = - 0.001 #* lt_kin ** 4 * ln_hd ** 1

    D02 = + 1. #* lt_kin ** 0 * ln_hd ** 2
    D12 = - 1. #* lt_kin ** 1 * ln_hd ** 2
    D22 = + 1. #* lt_kin ** 2 * ln_hd ** 2
    D32 = - 0.1 #* lt_kin ** 3 * ln_hd ** 2
    D42 = + 0.01 #* lt_kin ** 4 * ln_hd ** 2

    D03 = - 0.01 #* lt_kin ** 0 * ln_hd ** 3
    D13 = + 0.1 #* lt_kin ** 1 * ln_hd ** 3
    D23 = - 0.01  # * lt_kin ** 2 * ln_hd ** 3
    D33 = + 0.01  # * lt_kin ** 3 * ln_hd ** 3
    D43 = - 0.001  # * lt_kin ** 4 * ln_hd ** 3

    D04 = + 0.001  # * lt_kin ** 0 * ln_hd ** 4
    D14 = - 0.001  # * lt_kin ** 1 * ln_hd ** 4
    D24 = + 0.001 #* lt_kin ** 2 * ln_hd ** 4
    D34 = - 0.0001 #* lt_kin ** 3 * ln_hd ** 4
    D44 = + 0.00001 #* lt_kin ** 4 * ln_hd ** 4)


    D = D00,  D10, D20, D30, D40, \
        D01,  D11, D21, D31, D41, \
        D02,  D12, D22, D32, D42, \
        D03,  D13, D23, D33, D43, \
        D04,  D14, D24, D34, D44


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(x_fit, y_fit, data_to_fit)

    ax.set_xlabel('Temperature [K]')
    ax.set_ylabel('Density')
    ax.set_zlabel('cooling function')

    plt.show()

    #lambda_flatten = [item for sublist in data_to_fit for item in sublist]

    popt, pcov = \
    curve_fit(func, (np.log10(x_fit), np.log10(y_fit)),
                     np.log10(data_to_fit), D)

    return popt, pcov






