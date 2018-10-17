# -*- coding: utf-8 -*-
"""
Calculate equilibrium population of species and the cooling rate using
the Lipovka data for a certain kinetic temperatue and gas density

this script is mainly used for testing and development purposes
"""
from astropy import units as u

from frigus.cooling_function.fits import fit_lipovka
from frigus.population import cooling_rate_at_steady_state
from frigus.readers.dataset import DataLoader

species_data = DataLoader().load('HD_lipovka')

# density of the colliding species, in m^3, kinetic and radiation temperatures
nc_h = 1.e12 * u.meter ** -3
t_kin = 1200.0 * u.Kelvin
t_rad = 0.0 * u.Kelvin


cooling_rate = cooling_rate_at_steady_state(
    species_data,
    t_kin,
    t_rad,
    nc_h
).to(u.Joule / u.second)

cooling_rate_lipovka = fit_lipovka(t_kin, nc_h).to(u.Joule / u.second)

print(
    f'cooling rate\n'
    f' computed with frigus                           = {cooling_rate}\n'
    f' computed from the lipovka cooling function fit = {cooling_rate_lipovka}'
)

print('done')
