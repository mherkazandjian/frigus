# -*- coding: utf-8 -*-

#    metadata.py is part of Frigus.

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

"""Project metadata

Information describing the project.
"""


# The package name, which is also the "UNIX name" for the project.
package = 'frigus'
project = "level population and cooling function computation."
project_no_spaces = project.replace(' ', '')
version = '4.0.1'
description = (
    'Compute equilibrium and time-dependent popluation level densities\n'
    'for multi-level species. The level populations can be used to compute\n'
    'cooling function. Frigus allow for computing cooling functions\n'
    'parameterized by gas density, kinetic temperature and radiation\n'
    'temperatures.'
)
authors = ['Mher Kazandjian', 'Carla Maria Coppola']
authors_string = ', '.join(authors)
emails = ['mherkazandjian@gmail.com', 'carla.coppola@uniba.it']
license = 'GLP-3.0'
copyright = '2018 ' + authors_string
url = 'https://github.com/mherkazandjian/frigus'
download_url = 'https://github.com/mherkazandjian/frigus/archive/master.zip'
