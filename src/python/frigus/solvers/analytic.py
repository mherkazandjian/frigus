# -*- coding: utf-8 -*-

#    analytic.py is part of Frigus.

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

from numpy import exp, abs
from astropy.constants import k_B as kb


def population_density_ratio_two_level(g,
                                       energy_levels,
                                       k_10,
                                       a_10,
                                       n_c,
                                       t_kin,
                                       t_rad):
    """
    Calculate the equilibrium population density ratio of a two level system

    The provided parameters should all the compatible dimensions wise. i.e
    no checks are done if the input arguments

    At equilibrium, the ratio :math:`n_1 / n_0`:

   .. math::

        \\frac{n_1}{n_0} = \frac{B_{01} + K_{01}n_c}{A_{10} + B_{10} + K_{10}n_c} = \frac{M_{01}}{M_{10}}

    :param ndarray g: the degeneracies of the energy levels
    :param astropy.units.quantity.Quantity energy_levels: The energy levels
    :param astropy.units.quantity.Quantity k_10: The upper to lower
     collisional coefficient
    :param astropy.units.quantity.Quantity a_10: The upper to lower spontaneous
     emission rate
    :param astropy.units.quantity.Quantity n_c: The number density of the
     colliding species
    :param astropy.units.quantity.Quantity t_kin: The kintic temperature
    :param astropy.units.quantity.Quantity t_rad: The radiation temperature
    :return: float: The ratio of the upper to lower population density
    """
    g_0, g_1 = g
    delta_e = abs(energy_levels[1] - energy_levels[0])

    k_01 = (g_1 / g_0) * k_10 * exp(-delta_e / (kb * t_kin))

    f_10 = 1.0 / (exp(delta_e / (kb * t_rad)) - 1.0)

    b_10 = f_10 * a_10
    b_01 = (g_1 / g_0) * f_10 * a_10

    n_1_over_n_0 = (b_01 + k_01*n_c) / (a_10 + b_10 + k_10 * n_c)

    return n_1_over_n_0


def population_density_ratio_three_level(g,
                                         energy_levels,
                                         k_10,
                                         k_20,
                                         k_21,
                                         a_10,
                                         a_20,
                                         a_21,
                                         n_c,
                                         t_kin,
                                         t_rad):
    """
    Calculate the equilibrium population density ratio of a three level system

    The provided parameters should all the compatible dimensions wise. i.e
    no checks are done if the input arguments

    At equilibrium, the ratio :math:`n_1 / n_0`:

   .. math::

        \\frac{n_1}{n_0} = \frac{B_{01} + K_{01}n_c}{A_{10} + B_{10} + K_{10}n_c} = \frac{M_{01}}{M_{10}}

    :param ndarray g: the degeneracies of the energy levels
    :param astropy.units.quantity.Quantity energy_levels: The energy levels
    :param astropy.units.quantity.Quantity k_10: The upper to lower
     collisional coefficient from level 1 to level 0
    :param astropy.units.quantity.Quantity k_20: The upper to lower
     collisional coefficient from level 2 to level 0
    :param astropy.units.quantity.Quantity k_21: The upper to lower
     collisional coefficient from level 2 to level 1
    :param astropy.units.quantity.Quantity a_10: The upper to lower spontaneous
     emission rate from level 1 to level 0
    :param astropy.units.quantity.Quantity a_20: The upper to lower spontaneous
     emission rate from level 2 to level 0
    :param astropy.units.quantity.Quantity a_21: The upper to lower spontaneous
     emission rate from level 2 to level 1
    :param astropy.units.quantity.Quantity n_c: The number density of the
     colliding species
    :param astropy.units.quantity.Quantity t_kin: The kintic temperature
    :param astropy.units.quantity.Quantity t_rad: The radiation temperature
    :return: float: The ratio of the upper to lower population density
    """
    g_0, g_1, g_2 = g

    k_01 = (g_1 / g_0) * k_10 * exp(- abs(energy_levels[1] - energy_levels[0]) / (kb * t_kin))
    k_02 = (g_2 / g_0) * k_20 * exp(- abs(energy_levels[2] - energy_levels[0]) / (kb * t_kin))
    k_12 = (g_2 / g_1) * k_21 * exp(- abs(energy_levels[2] - energy_levels[1]) / (kb * t_kin))

    f_10 = 1.0 / (exp(abs(energy_levels[1] - energy_levels[0]) / (kb * t_rad)) - 1.0)
    f_20 = 1.0 / (exp(abs(energy_levels[2] - energy_levels[0]) / (kb * t_rad)) - 1.0)
    f_21 = 1.0 / (exp(abs(energy_levels[2] - energy_levels[1]) / (kb * t_rad)) - 1.0)

    b_10 = f_10 * a_10
    b_20 = f_20 * a_20
    b_21 = f_21 * a_21
    b_01 = (g_1 / g_0) * f_10 * a_10
    b_02 = (g_2 / g_0) * f_20 * a_20
    b_12 = (g_2 / g_1) * f_21 * a_21

    r_10 = k_10 * n_c + a_10 + b_10
    r_20 = k_20 * n_c + a_20 + b_20
    r_21 = k_21 * n_c + a_21 + b_21
    r_01 = k_01*n_c + b_01
    r_02 = k_02*n_c + b_02
    r_12 = k_12*n_c + b_12

    n_1_over_n_0 = ((r_01*r_20 + r_01*r_21 + r_21*r_02) /
                    (r_10*r_20 + r_10*r_21 + r_12*r_20))
    n_2_over_n_0 = ((r_02*r_10 + r_02*r_12 + r_12*r_01) /
                    (r_10*r_20 + r_10*r_21 + r_12*r_20))

    return n_1_over_n_0, n_2_over_n_0


def population_denisty_ratio_two_level_no_radiation(
        species_data,
        t_kin,
        n_collider):
    """
    Compute the analytic population density ration for a two level system

    The system is assumed to have no radiation field

    :param DataSetBase species_data: The dataset of a two level system
    :param t_kin: The kinetic temperature
    :param n_collider: The density of the colliding species
    :return: fractional population densities of the levels
    """
    a_10 = species_data.a_matrix[1, 0]
    k_10 = species_data.k_dex_matrix_interpolator(t_kin)[1, 0]
    n_c = n_collider
    g_0, g_1 = species_data.energy_levels.data['g']
    en_0, en_1 = species_data.energy_levels.data['E']

    k_01 = (g_1/g_0)*k_10*exp(-fabs(en_1 - en_0) / (kb * t_kin))

    r = n_c * k_01 / (n_c * k_10 + a_10)
    x_0 = 1.0 / (1.0 + r)
    x_1 = r / (1.0 + r)

    return x_0, x_1
