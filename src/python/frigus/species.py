# -*- coding: utf-8 -*-

#    dataset.py is part of Frigus.

#    Frigus: software to compure the energy exchange in a multi-level system
#    Copyright (C) 2016-2018 Mher V. Kazandjian and Carla Maria Coppola

#    Frigus is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#    Frigus is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Frigus.  If not, see <http://www.gnu.org/licenses/>.

"""
Module that implements species related classes e.g. energy levels
"""
import numpy
from astropy.table import QTable
from frigus import utils


class EnergyLevelsSpeciesBase(object):
    """
    Container class that holds energy levels data.
    """
    def __init__(self, n_levels=None, energy_unit=None, level_names=None,
                 *args, **kwargs):
        """
        Constructor
        """
        self.data = None
        """
        the astropy table that holds the levels data
        """

    def set_labels(self):
        """
        Set the field self.data['label']. This method modifies
        self.data['label'].
        """
        pass


class EnergyLevelsOnDegreeOfFreedom(EnergyLevelsSpeciesBase):
    """
    Energy levels of a one degree of freedom systems
    """
    def __init__(self,
                 n_levels=None,
                 energy_unit=None,
                 level_name='j',
                 *args, **kwargs):
        """
        Constructor
        """
        super(EnergyLevelsOnDegreeOfFreedom, self).__init__(*args, **kwargs)

        table_type = {
            level_name: numpy.zeros(n_levels, 'i4'),
            'g': numpy.zeros(n_levels, 'i4'),
            'label': numpy.zeros(n_levels, 'i4'),
            'E': numpy.zeros(n_levels, 'f8') * energy_unit
        }
        self.data = QTable(dict(**table_type))
        """the astropy table that holds the levels data"""


class EnergyLevelsMolecule(EnergyLevelsSpeciesBase):
    """
    Container class that holds energy levels data. This supports two quantum
    numbers of the molecule, vibrations and rotations.

    .. todo:: rename this class to EnergyLevelsTwoDegreesOfFreedom
    """
    def __init__(self, n_levels, energy_unit=None, *args, **kwargs):
        """
        Constructor
        """
        super(EnergyLevelsMolecule, self).__init__(n_levels, *args, **kwargs)

        table_type = {
            'v': numpy.zeros(n_levels, 'i4'),
            'j': numpy.zeros(n_levels, 'i4'),
            'g': numpy.zeros(n_levels, 'i4'),
            'label': numpy.zeros(n_levels, 'i4'),
            'E': numpy.zeros(n_levels, 'f8') * energy_unit
        }
        self.data = QTable(dict(**table_type))
        """the astropy table that holds the levels data"""

        self.v_max_allowed = None
        """The maximum value of v that is allowed. This could be higher or
        lower than the values in self.data['v']"""

        self.j_max_allowed = None
        """The maximum value of j that is allowed. This could be higher or
        lower than the values in self.data['j']"""

    def set_labels(self, v_max=None):
        """
        set the field self.data['label']. This method modifies
        self.data['label'] and self.v_max_allowed

        :param v_max: The maximum value of v to be used in computing and
         setting the labels.
        """
        if v_max is not None:
            self.v_max_allowed = v_max
        self.data['label'] = utils.linear_2d_index(self.data['v'],
                                                   self.data['j'],
                                                   n_i=v_max)
