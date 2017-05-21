# -*- coding: utf-8 -*-
"""
calculate equilibrium population of species and the cooling function.
"""
from __future__ import print_function
import pylab
import numpy

from astropy import units as u
import matplotlib.pyplot as plt

from population import population_density_at_steady_state
import utils

from readers import DataLoader

species_data = DataLoader().load('H2_lique')
# species_data = DataLoader().load('HD_lipovka')
def ode_rhs(t, y, args):
    """defining the function which will be the rhs of df/dt"""
    return dot(full, y)

# evolving with respect to time using scipy integrator
if True:
    r = ode(ode_rhs, jac = None).set_integrator('vode',
    # r = ode(ode_rhs, jac = None).set_integrator('dopri',
                                                method='bdf',
                                                #with_jacobian = False,
                                                rtol = 1e-3)
    r.set_initial_value(f0, t0).set_f_params(1.0)

    i = 0
    # while r.successful() and r.t < tf:
    while r.t < tf:
        r.integrate(r.t + dt)
        t.append(r.t)
        ft_0.append(r.y[lPlot[0]])
        ft_1.append(r.y[lPlot[1]])
        ft_2.append(r.y[lPlot[2]])
        ft_3.append(r.y[lPlot[3]])
        ft_4.append(r.y[lPlot[4]])
        ft_5.append(r.y[lPlot[5]])

        if i % 100 == 0:
            print 'i = %d t = %e' % (i, r.t), 1.0 - numpy.sum(r.y), 1.0 - r.y[0]/f[0]
            print 'stepNum = %d t = %e, current dt = %e, 1 - sum(pop_dens) = %e' % (i, r.t, dt, 1.0 - numpy.sum(r.y))

        i += 1

# solving using my implementation of the bulrische stoer integrator
if False:
    rel_tol = 1e-3
    dt0_bs = 1.0

    # initial conditions
    state = mylib.numerics.ode.State(f0.size)
    state.t = 0.0
    state.y[:] = f0

    # defining the function which will be the rhs of df/dt
    def ode_rhs(state, state_der, parms=None):

        state_der.y[:] = dot(full, state.y)
        return True

    solver = mylib.numerics.ode.BS(log_dir='.', dx0=dt0_bs, reltol=rel_tol,
                                   state0=state, verbose=False)
    solver.set_derivs_func(ode_rhs)

    stepNum = 0
    while solver.state.x <= tf:

        t.append(solver.state.x)
        ft_0.append(solver.state.y[0])
        ft_1.append(solver.state.y[1])
        ft_2.append(solver.state.y[2])
        ft_3.append(solver.state.y[3])
        ft_4.append(solver.state.y[4])
        ft_5.append(solver.state.y[5])

        solver.advance_one_step()

        if stepNum % 100 == 0:
            print 'stepNum = %d t = %e, current dt = %e, 1 - sum(pop_dens) = %e' %\
                  (stepNum, solver.state.x, solver.dx,
                   1.0 - numpy.sum(solver.state.y))

        stepNum += 1

print('done')