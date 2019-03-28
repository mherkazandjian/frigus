# -*- coding: utf-8 -*-

#    fits.py is part of Frigus.

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
"""
module that implements cooling functions or fits of cooling functions
"""
import numpy
from astropy import units as u


def fit_glover(t_kin):
    """
    Compute the cooling rate of the Glover Fit of the cooling rate of H2

    The cooling rate is computed in units of erg / cm^3 / s^-1

    reference: https://goo.gl/htEPuJ

    :param ndarray|float t_kin: The input kinetic temperature (in K)
    :return: ndarray|float the cooling rate
    """
    if 100.0 <= t_kin <= 1000.0:
        retval = 10**(-24.311209
                      + 3.5692468*numpy.log10(t_kin/1000.)
                      - 11.332860*(numpy.log10(t_kin/1000.))**2
                      - 27.850082*(numpy.log10(t_kin/1000.))**3
                      - 21.328264*(numpy.log10(t_kin/1000.))**4
                      - 4.2519023*(numpy.log10(t_kin/1000.))**5)
    elif 1000.0 < t_kin <= 6000.0:
        retval = 10**(-24.311209
                      + 4.6450521*numpy.log10(t_kin/1000.)
                      - 3.7209846*numpy.log10((t_kin/1000.))**2
                      + 5.9369081*numpy.log10((t_kin/1000.))**3
                      - 5.5108047*numpy.log10((t_kin/1000.))**4
                      + 1.5538288*numpy.log10((t_kin/1000.))**5)
    else:
        msg = "inpute kinetic temperature {} is out of bounds".format(t_kin)
        raise ValueError(msg)

    return retval * u.erg * u.s**-1 * u.cm**3


def fit_lipovka(t_kin, n_hd):
    """
    Compute the cooling rate of  HD as a function of temperature and density

    The cooling function is computed in units of erg / s

    reference: HD revised cooling function - Lipovka, Nunez, Avila Reese 2005

    :param u.quantity.Quantity t_kin: The temperature range as an astropy
     quantity (in K)
    :param u.quantity.Quantity n_hd: The density of HD as an astropy quantity
    :return: u.quantity.Quantity: The cooling function in erg / s for the
     inpute T values

    .. todo:: add an example script for plotting this cooling function
    """
    #
    # make sure the input temperatures and density values are within the
    # validity range of the fit and are of the correct types
    #
    assert isinstance(t_kin, u.quantity.Quantity)
    assert isinstance(n_hd, u.quantity.Quantity)

    assert t_kin.min() >= 100.0 * u.K
    assert t_kin.max() <= 2000.0 * u.K
    assert n_hd.min() >= 1e6 * u.m**-3
    assert n_hd.max() <= 1e14 * u.m**-3

    lt_kin = numpy.log10(t_kin.value)
    ln_hd = numpy.log10(n_hd.cgs.value)

    retval = 10.**(- 42.57688 * lt_kin ** 0 * ln_hd ** 0
                   + 0.92433 * lt_kin ** 0 * ln_hd ** 1
                   + 0.54962 * lt_kin ** 0 * ln_hd ** 2
                   - 0.07676 * lt_kin ** 0 * ln_hd ** 3
                   + 0.00275 * lt_kin ** 0 * ln_hd ** 4

                   + 21.93385 * lt_kin ** 1 * ln_hd ** 0
                   + 0.77952 * lt_kin ** 1 * ln_hd ** 1
                   - 1.06447 * lt_kin ** 1 * ln_hd ** 2
                   + 0.11864 * lt_kin ** 1 * ln_hd ** 3
                   - 0.00366 * lt_kin ** 1 * ln_hd ** 4

                   - 10.19097 * lt_kin ** 2 * ln_hd ** 0
                   - 0.54263 * lt_kin ** 2 * ln_hd ** 1
                   + 0.62343 * lt_kin ** 2 * ln_hd ** 2
                   - 0.07366 * lt_kin ** 2 * ln_hd ** 3
                   + 0.002514 * lt_kin ** 2 * ln_hd ** 4

                   + 2.19906 * lt_kin ** 3 * ln_hd ** 0
                   + 0.11711 * lt_kin ** 3 * ln_hd ** 1
                   - 0.13768 * lt_kin ** 3 * ln_hd ** 2
                   + 0.01759 * lt_kin ** 3 * ln_hd ** 3
                   - 0.000666317 * lt_kin ** 3 * ln_hd ** 4

                   - 0.17334 * lt_kin ** 4 * ln_hd ** 0
                   - 0.00835 * lt_kin ** 4 * ln_hd ** 1
                   + 0.0106 * lt_kin ** 4 * ln_hd ** 2
                   - 0.001482 * lt_kin ** 4 * ln_hd ** 3
                   + 0.000061926 * lt_kin ** 4 * ln_hd ** 4)

    return retval * u.erg * u.s**-1


def fit_lipovka_low_density(t_kin):
    """
    Compute the fit of the cooling rate of HD at low kinetic temperature(s)

    The cooling function is computed for low kinetic temperatures in units of
    of erg / s for low temperature

    reference: HD revised cooling function - Lipovka, Nunez, Avila Reese 2005

    :param ndarray|float t_kin: The kinetic temperature at which the cooling
     function will be computedd in K
    :return: ndarray|float: the cooling function in units of erg / s
    """
    assert isinstance(t_kin, u.quantity.Quantity)

    assert t_kin.min() >= 100.0 * u.K
    assert t_kin.max() <= 20000.0 * u.K

    lt_kin = numpy.log10(t_kin.value)

    retval = 10.0**(- 42.45906
                    + 21.90083 * lt_kin
                    - 10.1954 * lt_kin**2
                    + 2.19788 * lt_kin**3
                    - 0.17286 * lt_kin**4)

    return retval * u.erg * u.s**-1


def fit_coppola(t_kin, n_hd):
    """
    Compute the cooling rate of H2 as a function of temperature and density

    The cooling function is computed in units of erg / s

    reference: Updated temperature and density dependent cooling function for
               H$_2$ colliding with atomic hydrogen. -
               Coppola, Lique, Mazzia & Kazandjian 2018

    :param u.quantity.Quantity t_kin: The temperature range as an astropy
     quantity (in K)
    :param u.quantity.Quantity n_hd: The density of HD as an astropy quantity
    :return: u.quantity.Quantity: The cooling function in erg / s for the
     input T values

    """
    #
    # make sure the input temperatures and density values are within the
    # validity range of the fit and are of the correct types
    #
    assert isinstance(t_kin, u.quantity.Quantity)
    assert isinstance(n_hd, u.quantity.Quantity)

    assert t_kin.min() >= 100.0 * u.K
    assert t_kin.max() <= 5000.0 * u.K
    assert n_hd.min() >= 1e2 * u.m**-3
    assert n_hd.max() <= 1e14 * u.m**-3

    lt_kin = numpy.log10(t_kin.value)
    ln_hd = numpy.log10(n_hd.cgs.value)

    # retval = 10.**(- 2.82483125e2 * lt_kin ** 0 * ln_hd ** 0
    #                - 1.34604375e2 * lt_kin ** 0 * ln_hd ** 1
    #                + 1.07432776e2 * lt_kin ** 0 * ln_hd ** 2
    #                - 2.04446787e1 * lt_kin ** 0 * ln_hd ** 3
    #                + 1.17030794e0 * lt_kin ** 0 * ln_hd ** 4
    #
    #                + 3.92355880e2 * lt_kin ** 1 * ln_hd ** 0
    #                + 2.15639291e2 * lt_kin ** 1 * ln_hd ** 1
    #                - 1.68827762e2 * lt_kin ** 1 * ln_hd ** 2
    #                + 3.18279345e1 * lt_kin ** 1 * ln_hd ** 3
    #                - 1.80956442e0 * lt_kin ** 1 * ln_hd ** 4
    #
    #                - 2.26543263e2 * lt_kin ** 2 * ln_hd ** 0
    #                - 1.26554822e2 * lt_kin ** 2 * ln_hd ** 1
    #                + 9.80765328e1 * lt_kin ** 2 * ln_hd ** 2
    #                - 1.83570655e1 * lt_kin ** 2 * ln_hd ** 3
    #                + 1.03763936e0 * lt_kin ** 2 * ln_hd ** 4
    #
    #                + 5.80077743e1 * lt_kin ** 3 * ln_hd ** 0
    #                + 3.25577420e1 * lt_kin ** 3 * ln_hd ** 1
    #                - 2.50129961e1 * lt_kin ** 3 * ln_hd ** 2
    #                + 4.65210999e0 * lt_kin ** 3 * ln_hd ** 3
    #                - 2.61574544e-1 * lt_kin ** 3 * ln_hd ** 4
    #
    #                - 5.50713938e0 * lt_kin ** 4 * ln_hd ** 0
    #                - 3.10446347e0 * lt_kin ** 4 * ln_hd ** 1
    #                + 2.36722261e0 * lt_kin ** 4 * ln_hd ** 2
    #                - 4.37732460e-1 * lt_kin ** 4 * ln_hd ** 3
    #                + 2.44900752e-2 * lt_kin ** 4 * ln_hd ** 4)

    # retval = 10.**(  -9.30186394e+01 * lt_kin ** 0 * ln_hd ** 0
    #                  -4.59042326e+00* lt_kin ** 0 * ln_hd ** 1
    #                  -2.39644612e+00* lt_kin ** 0 * ln_hd ** 2
    #                  +5.33417509e-01* lt_kin ** 0 * ln_hd ** 3
    #                  -2.37275024e-02* lt_kin ** 0 * ln_hd ** 4
    # 
    #                  +9.53790916e+01 * lt_kin ** 1 * ln_hd ** 0
    #                  +1.21991916e+01* lt_kin ** 1 * ln_hd ** 1
    #                  +2.84793324e-01* lt_kin ** 1 * ln_hd ** 2
    #                  -1.79478025e-01* lt_kin ** 1 * ln_hd ** 3
    #                  -5.06173402e-04* lt_kin ** 1 * ln_hd ** 4
    # 
    #                  -5.29347560e+01* lt_kin ** 2 * ln_hd ** 0
    #                  -8.19882323e+00* lt_kin ** 2 * ln_hd ** 1
    #                  +1.18055362e+00* lt_kin ** 2 * ln_hd ** 2
    #                  -1.88445938e-01* lt_kin ** 2 * ln_hd ** 3
    #                  +1.82380851e-02* lt_kin ** 2 * ln_hd ** 4
    # 
    #                  +1.32856875e+01* lt_kin ** 3 * ln_hd ** 0
    #                  +2.20881654e+00* lt_kin ** 3 * ln_hd ** 1
    #                  -5.29631603e-01* lt_kin ** 3 * ln_hd ** 2
    #                  +1.03933459e-01* lt_kin ** 3 * ln_hd ** 3
    #                  -8.27289914e-03* lt_kin ** 3 * ln_hd ** 4
    # 
    #                  -1.22881695e+00* lt_kin ** 4 * ln_hd ** 0
    #                  -2.10210825e-01* lt_kin ** 4 * ln_hd ** 1
    #                  +6.52317292e-02* lt_kin ** 4 * ln_hd ** 2
    #                  -1.39994365e-02* lt_kin ** 4 * ln_hd ** 3
    #                  +1.06496239e-03* lt_kin ** 4 * ln_hd ** 4)

    retval = 10.**(  -1.07761178e+02 * lt_kin ** 0 * ln_hd ** 0
                     +8.50704285e+00* lt_kin ** 0 * ln_hd ** 1
                     -3.08850636e-01* lt_kin ** 0 * ln_hd ** 2
                     -3.97357274e-01* lt_kin ** 0 * ln_hd ** 3
                     +3.65589231e-02* lt_kin ** 0 * ln_hd ** 4

                     +1.17901741e+02 * lt_kin ** 1 * ln_hd ** 0
                     -1.09489742e+01* lt_kin ** 1 * ln_hd ** 1
                     +2.51708582e-01* lt_kin ** 1 * ln_hd ** 2
                     +5.47452438e-01* lt_kin ** 1 * ln_hd ** 3
                     -4.89055533e-02* lt_kin ** 1 * ln_hd ** 4

                     -6.61721218e+01* lt_kin ** 2 * ln_hd ** 0
                     +5.95232129e+00* lt_kin ** 2 * ln_hd ** 1
                     -1.10639748e-01* lt_kin ** 2 * ln_hd ** 2
                     -2.92423859e-01* lt_kin ** 2 * ln_hd ** 3
                     +2.57375696e-02* lt_kin ** 2 * ln_hd ** 4

                     +1.67152116e+01* lt_kin ** 3 * ln_hd ** 0
                     -1.43513388e+00* lt_kin ** 3 * ln_hd ** 1
                     +2.84953800e-02* lt_kin ** 3 * ln_hd ** 2
                     +7.02040444e-02* lt_kin ** 3 * ln_hd ** 3
                     -6.18478233e-03* lt_kin ** 3 * ln_hd ** 4

                     -1.55475014e+00* lt_kin ** 4 * ln_hd ** 0
                     +1.29639563e-01* lt_kin ** 4 * ln_hd ** 1
                     -3.13943372e-03* lt_kin ** 4 * ln_hd ** 2
                     -6.36261310e-03* lt_kin ** 4 * ln_hd ** 3
                     +5.65797161e-04* lt_kin ** 4 * ln_hd ** 4)


    return retval * u.erg * u.s**-1


def fit_flower(density, t_kin):
    """
    Compute the cooling function by Flower et at for HD colliding with H
    as reported in the subroutine
    reference: http://massey.dur.ac.uk/drf/cooling_by_HD/W_HD.f

    The cooling rate is computed in units of erg / s^-1

    :param u.quantity.Quantity t_kin: temperature range as an astropy
           quantity (in K)
    :param u.quantity.Quantity density: density as an astropy quantity (in cm-3)
    :return: ndarray|float the cooling rate
    """

    assert isinstance(t_kin, u.quantity.Quantity)
    assert isinstance(density, u.quantity.Quantity)

    aa = -26.2982
    bb = -0.215807
    omeg = 2.03657
    phi = 4.63258
    c1 = 0.283978
    c2 = -1.27333
    d1 = -2.08189
    d2 = 4.66288

    y = numpy.log10(density.cgs.value)
    x = numpy.log10(t_kin.si.value)

    w = 0.5 * y + aa * x**bb - \
        numpy.sqrt(
            0.25 * y*y + (c1*y+c2) * numpy.sin(omeg*x+phi) + (d1*y+d2)
        )

    retval = 10.0**w

    return retval * u.erg * u.s**-1


def fit_associative_detachment(t_kin):
    """
    Compute the reaction rate by Coppola et al 2007 for H-/H associative
    detachment

    The reaction rate is computed in units of m**3 s**-1

    :param u.quantity.Quantity t_kin: temperature range as an astropy
           quantity (in K)
    :return: u.quantity |float the reaction rate
    """

    # nc_hm = nc_h * 1e-10
    rr_coppola2007 = 3.67e-5 \
                 * (t_kin.value ** (-2.28)) \
                 * numpy.exp(-47172. / t_kin.value)

    # m_h_hm = rr_coppola2007 * nc_h * nc_hm * numpy.ones(55) /55.

    return rr_coppola2007 * u.m ** 3 * u.s ** -1
