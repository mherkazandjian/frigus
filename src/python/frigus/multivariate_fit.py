import numpy as np
from scipy.optimize import curve_fit

def func(X, D):
    T, nh = X

    D00,  D10, D20, D30, D40, \
    D01,  D11, D21, D31, D41, \
    D02,  D12, D22, D32, D42, \
    D03,  D13, D23, D33, D43, \
    D04,  D14, D24, D34, D44 = D

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
    lin12 = D31 * (T**2) * nh
    lin13 = D41 * (T**3) * nh
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