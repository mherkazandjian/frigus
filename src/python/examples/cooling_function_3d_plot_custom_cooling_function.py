"""
compute the cooling function as a function of kinetic temperature and gas
density and plot them in 3D
"""
from __future__ import print_function
import numpy as np
from astropy import units as u
from frigus.readers.dataset import DataLoader
from frigus.cooling_function.grid import CoolingFunctionGrid
from frigus.population import compute_transition_rate_matrix, cooling_rate
from frigus.solvers.linear import solve_equilibrium


class CoolingFunctionGridMultiSpecies(CoolingFunctionGrid):
    def compute(self):
        self._compute_mesh()

        # ----------------------------------------------------------------------
        species_data_h2_h = DataLoader().load('H2_lique')

        # load the data of He + H2 w/wo  radiative data
        species_data_h2_he = DataLoader().load('HeH2')
        species_data_h2_he_zero_a = DataLoader().load('HeH2')
        species_data_h2_he_zero_a.a_matrix *= 0.0
        species_data_h2_he_zero_a.raw_data.a *= 0.0
        species_data_h2_he_zero_a.raw_data.a_info_nnz = list(
            species_data_h2_he_zero_a.raw_data.a_info_nnz
        )
        species_data_h2_he_zero_a.raw_data.a_info_nnz[4] *= 0.0

        # load the data of H+ + H2 w/wo  radiative data
        species_data_h2_hp = DataLoader().load('HpH2')
        species_data_h2_hp_zero_a = DataLoader().load('HpH2')
        species_data_h2_hp_zero_a.a_matrix *= 0.0
        species_data_h2_hp_zero_a.raw_data.a *= 0.0
        species_data_h2_hp_zero_a.raw_data.a_info_nnz = list(
            species_data_h2_hp_zero_a.raw_data.a_info_nnz
        )
        species_data_h2_hp_zero_a.raw_data.a_info_nnz[4] *= 0.0
        # ----------------------------------------------------------------------

        cooling_rate_grid = np.zeros_like(self.n_grid.value).flatten()
        for i, (n, t_kin, t_rad) in enumerate(
                zip(self.n_grid.flat,
                    self.t_kin_grid.flat,
                    self.t_rad_grid.flat)):

            # -------------------------------------
            m_h2_h = compute_transition_rate_matrix(species_data_h2_h,
                                                    t_kin,
                                                    t_rad,
                                                    n)

            m_h2_he_without_a = compute_transition_rate_matrix(
                species_data_h2_he_zero_a,
                t_kin,
                t_rad,
                n * 1.e-1)

            m_h2_hp_without_a = compute_transition_rate_matrix(
                species_data_h2_hp_zero_a,
                t_kin,
                t_rad,
                n* 2.e-4)

            m = m_h2_h + m_h2_he_without_a + m_h2_hp_without_a

            x_equilibrium = solve_equilibrium(m.si.value)

            sca = cooling_rate(
                x_equilibrium,
                species_data_h2_h.energy_levels,
                species_data_h2_h.a_matrix
            )

            cooling_rate_grid[i] = sca.cgs.value

        cooling_rate_grid = cooling_rate_grid.reshape(self.n_grid.shape)
        self.cooling_function = cooling_rate_grid
        return cooling_rate_grid


grid = CoolingFunctionGridMultiSpecies()

grid.set_species(DataLoader().load('HD_lipovka'))

grid.set_density(np.logspace(6, 14, 35) * u.m ** -3)
grid.set_t_kin(np.logspace(2, 3.2, 35) * u.Kelvin)
grid.set_t_rad(0.0 * u.Kelvin)

grid.plot_3d(x='n', y='t_kin', show=True)

print('done')
