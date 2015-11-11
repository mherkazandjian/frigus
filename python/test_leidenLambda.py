import sys, os
import matplotlib
matplotlib.use('Qt4Agg')
import pylab

import molData
import numpy
from numpy import float64, abs, zeros, ones, zeros, arange, exp, eye, dot
from numpy import linalg, vstack, array
from ismUtils import planckOccupation as ng
import mylib.numerics.ode

restore = False
# reading the whole database of line info of species from LAMBDA
lambdaPath = '/home/mher/ism/code/radex/Radex/data/home.strw.leidenuniv.nl/~moldata/datafiles'

# reads the whole database
# ------------------------
# reader = molData.reader(dirPath = lambdaPath)

# retrieves only one specie
reader = molData.reader(dirPath = lambdaPath, species = 'NH3')

Tkin_min = 15.0
Tkin_max = 300.0
TkinSamples = 40.0 
nH2 = [1e2, 1e4, 1e6]

Tcmb = 2.73  # temperature of background radiation (K)

# some constants
hPlank = 6.63e-27    # erg.s
cLight = 29979245800.0 # cm.s^-1
kBoltz = 1.38e-16    # erg.K^-1 
ev2erg = 1.602e-12   #erg

# selecting the one holding the info for p-NH3
pNH3 = reader.get_specie(specStr='NH3', inInfo='p-')

# converting the energy units to K 
for idx, level in enumerate(pNH3.levels):
    level['E'] *= hPlank*cLight/kBoltz # energies in K


def computeRateMatrix(pNH3, Tkin, nc):
    """compute the matrix of transition rates"""

    n = pNH3.nlevels

    levels = pNH3.levels
    transRad = pNH3.transRad
    transColl = pNH3.transColl['p-H2']['trans']
      
    ###########################################################
    # constructing the matrix 
    ###########################################################
    def fill_K_matrix():
        """fill the kij matrix from the lambda collisional transition
        database"""

        # n  = 5
        K = zeros( (n, n), dtype = float64)
        for trans in transColl:
            u  = trans['u']; l = trans['l']
            gu = levels[u]['g']; gl = levels[l]['g']

            # difference in the enrergy in K
            dE = abs(levels[u]['E'] - levels[l]['E'])
            
            K[u, l] = trans['rc'](Tkin)
            K[l, u] = (float64(gu) / float64(gl)) * K[u,l] * exp(-dE / Tkin)
        
        return K

    def fill_AP_matrix():
        """fill the (A prime)_ij matrix from the lambda radiative transitions
        database"""
        AP = zeros( (n, n), dtype = float64)
        
        for trans in transRad:
            u  = trans['u']; l = trans['l']
            dE = abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
            
            nu = dE*kBoltz / hPlank # freq in Hz
            
            AP[u, l] = (1.0 + ng(hPlank, nu, kBoltz, Tcmb))*trans['A']
    
        return AP
    
    def fill_ABS_matrix():
        """fill the Aij matrix for absorbtion transitions from the lambda
        radiative transitions database"""
        ABS = zeros( (n, n), dtype = float64)
        
        for trans in transRad:
            u  = trans['u']; l = trans['l']
            dE = abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
            gu, gl = float64(levels[u]['g']), float64(levels[l]['g'])
            
            nu = dE*kBoltz / hPlank # freq in Hz
            
            ABS[u, l] = (gu/gl)*ng(hPlank, nu, kBoltz, Tcmb)*trans['A']
    
        return ABS
    
    def fill_E_matrix():
        """This is a matrix full of ones (it is a utility matrix"""
        return ones((n,n))
    
    K = fill_K_matrix()
    AP = fill_AP_matrix()
    ABS = fill_ABS_matrix()
    E = fill_E_matrix()
    
    F = nc * K + AP + ABS.T
    diag = eye(n)*dot(F, E)[:, 1]
    offdiag = F.T
    
    full = -diag + offdiag

    return full


def solveEquilibrium(pNH3, full):
    """solve for the equlibrium population densities"""
    n = pNH3.nlevels

    # solving directly
    #replacing the first row with the conservation equation
    dndt = zeros((n,1))
    full[0,:] = 1.0
    dndt[0]   = 1.0
    
    A = full
    b = dndt
    #solving the system A.x = b
    #before solving, we will devide each row by the diagonal
    for i in arange(n):
        A[i,:] = A[i,:]/A[i,i]
    x = linalg.solve(A, b)
    
    #print x.T
    
    # the fractional population density
    f = x
    return f


def computeLuminosity(pNH3):
    """.. todo:: add doc"""
    n = pNH3.nlevels
    levels = pNH3.levels
    transRad = pNH3.transRad


    def fill_R_matrix():
        """fill the L_ij matrix from the lambda radiative transitions
        database"""
        R = zeros( (n, n), dtype = float64)
        
        for trans in transRad:
            u  = trans['u']; l = trans['l']
            dE = abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
            
            nu = dE*kBoltz / hPlank # freq in Hz
            
            R[u, l] = trans['A']*hPlank*nu
    
        return R

    R = fill_R_matrix()
    
    return R

####### solve for one set of parameters ######
Tkin = 30.0
nc = 1000.0

# solving using matrix inversion
full = computeRateMatrix(pNH3, Tkin, nc)
f = solveEquilibrium(pNH3, full.copy())

# solving using minimization
full2 = computeRateMatrix(pNH3, Tkin, nc)

# generating the constraint equation and appending it to the Matrix
cons  = ones(pNH3.nlevels)
full2 = vstack((full2, cons))

# generating the RHS
rhs   = zeros(pNH3.nlevels + 1)
rhs[-1] = 1

# solving
sol = linalg.lstsq(full2, rhs)
f2 = sol[0]

pylab.figure( figsize = (6,12))
pylab.subplot(211)
pylab.semilogy(f)
pylab.hold(True)
pylab.semilogy(f2, 'o')
pylab.xlabel('level')
pylab.ylabel('population density')

# solving by evolving with time
f0 = zeros(pNH3.nlevels)

# setting up the initial conditions to population densities
# that would be attained at LTE
if False:
    f0 = zeros(pNH3.nlevels)
    Z = float64(0.0)  # the partition function
    for i in arange(pNH3.nlevels):
        f0[i] = pNH3.levels[i]['g']*exp(- pNH3.levels[i]['E'] / Tkin)
        Z += f0[i]
    f0 /= Z # the initial fractional population densities
    
if False:
    f0 = f # using the eq sol as ICs

if True:
    # setting population levels manually
    f0[0] = 0.3
    f0[1] = 1e-2
    f0[2] = 1e-3
    f0[3] = 0.5
    f0[4] = 1e-5
    f0[5] = 1e-6
    
    # normalizing to 1
    f0 /= f0.sum() 
    
t0 = 0.0 # the initial time
dt = 2e4 # initial timestep in seconds
tf = 2e8 # final time

lPlot = array([0, 1, 2, 3, 4, 5])
t = []
ft_0 = []
ft_1 = []
ft_2 = []
ft_3 = []
ft_4 = []
ft_5 = []


def ode_rhs(t, y, args):
    """defining the function which will be the rhs of df/dt"""
    return dot(full, y)

'''
#evolving with respect to time
from scipy.integrate import ode
#r = ode(ode_rhs, jac = None).set_integrator('vode', 
r = ode(ode_rhs, jac = None).set_integrator('dopri', 
                                            method='bdf', 
                                            #with_jacobian = False,
                                            rtol = 1e-12)
r.set_initial_value(f0, t0).set_f_params(1.0)

i = 0
while r.successful() and r.t < tf:
    r.integrate(r.t+dt)
    t.append(r.t)
    ft_0.append(r.y[ lPlot[0] ])
    ft_1.append(r.y[ lPlot[1] ])
    ft_2.append(r.y[ lPlot[2] ])
    ft_3.append(r.y[ lPlot[3] ])
    ft_4.append(r.y[ lPlot[4] ])
    ft_5.append(r.y[ lPlot[5] ])

    if i % 100 == 0:
        print 'i = %d t = %e' % (i, r.t), 1.0 - numpy.sum(r.y), 1.0 - r.y[0]/f[0]
        print 'stepNum = %d t = %e, current dt = %e, 1 - sum(pop_dens) = %e' % (i, r.t, dt, 1.0 - numpy.sum(r.y))

    i+=1
'''

'''
####################################################################################################
            solving using my implementation of the bulrische stoer integrator
####################################################################################################
'''
rel_tol = 1e-3
n_steps = 1000
dt0_bs = 1.0

# intial conditions
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

####### finishing integrating using the BS method #########

t = array(t)

# pylab.figure(1)
pylab.subplot(212)
# plotting the actual curves with the equilib sols (dashes)
pylab.hold(True)
pylab.loglog(t, ft_0,'r')
pylab.loglog(t, ft_1,'g')
pylab.loglog(t, ft_2,'b')
# pylab.loglog(t, ft_3,'c')
# pylab.loglog(t, ft_4,'k')
# pylab.loglog(t, ft_5,'y')
pylab.loglog([dt,tf], [f2[lPlot[0]],f2[lPlot[0]]], '--r')
pylab.loglog([dt,tf], [f2[lPlot[1]],f2[lPlot[1]]], '--g')
pylab.loglog([dt,tf], [f2[lPlot[2]],f2[lPlot[2]]], '--b')
# pylab.loglog([dt,tf], [f2[lPlot[3]],f2[lPlot[3]]], '--c')
# pylab.loglog([dt,tf], [f2[lPlot[4]],f2[lPlot[4]]], '--k')
# pylab.loglog([dt,tf], [f2[lPlot[5]],f2[lPlot[5]]], '--y')
pylab.axis([dt, t.max()*1.1, 1e-13, 1])
pylab.xlabel('time (s)')
pylab.ylabel('population density')

"""
pylab.hold(True)
#plotting the relative difference between the final sol and the eq sol
pylab.loglog(t, numpy.fabs(1.0 - ft_0/f2[lPlot[0]]),'r')
pylab.loglog(t, numpy.fabs(1.0 - ft_1/f2[lPlot[1]]),'g')
pylab.loglog(t, numpy.fabs(1.0 - ft_2/f2[lPlot[2]]),'b')
pylab.axis([dt, tf, 1e-9, 1])
"""

pylab.ion()
pylab.show()

"""
####### plotting the line intensities per molecule as a function of Tkin#####
for nc in nH2:
    
    Tkins = numpy.linspace(Tkin_min, Tkin_max, TkinSamples)
    Llines = {'11':[], '22':[], '44':[]}
    
    for Tkin in Tkins:
        
        full = computeRateMatrix(pNH3, Tkin, nc)
        f    = solveEquilibrium(pNH3, full)
        R    = computeLuminosity(pNH3)
        
        u11 = 1; l11 = 0;    
        L11 = R[u11][l11]*f[u11]
        
        u22 = 3; l22 = 2;    
        L22 = R[u22][l22]*f[u22]
        
        u44 = 11; l44 = 10;    
        L44 = R[u44][l44]*f[u44]
    
        Llines['11'].append(L11)
        Llines['22'].append(L22)
        Llines['44'].append(L44)
        
    if nc == 100.0:
        colorStr = 'g'
    if nc == 10000.0:
        colorStr = 'r'
    if nc == 1000000.0:
        colorStr = 'b'
        
    pylab.plot(Tkins, Llines['11'], colorStr)
    pylab.hold(True)
    pylab.plot(Tkins, Llines['22'], colorStr+'--')
    pylab.plot(Tkins, Llines['44'], colorStr+'-.')
    
pylab.xscale('log')    
pylab.yscale('log')
pylab.axis([15,300,1e-25,1e-22])
pylab.show()
"""

"""
################ plotting the line ratios as a function of Tkin###############
for nc in nH2:
    
    Tkins = numpy.linspace(Tkin_min, Tkin_max, TkinSamples)
    Llines = {'11':[], '22':[], '44':[]}
    
    for Tkin in Tkins:
        
        full = computeRateMatrix(pNH3, Tkin, nc)
        f    = solveEquilibrium(pNH3, full)
        R    = computeLuminosity(pNH3)
        
        u11 = 1; l11 = 0;    
        L11 = R[u11][l11]*f[u11]
        
        u22 = 3; l22 = 2;    
        L22 = R[u22][l22]*f[u22]
        
        u44 = 11; l44 = 10;    
        L44 = R[u44][l44]*f[u44]
    
        Llines['11'].append(L11)
        Llines['22'].append(L22)
        Llines['44'].append(L44)
        
    if nc == 100.0:
        colorStr = 'g'
    if nc == 10000.0:
        colorStr = 'r'
    if nc == 1000000.0:
        colorStr = 'b'
        
    pylab.plot(Tkins, array(Llines['22'])/array(Llines['11']), colorStr)
    pylab.hold(True)
    pylab.plot(Tkins, array(Llines['44'])/array(Llines['22']), colorStr+'--')
    pylab.plot(Tkins, array(Llines['44'])/array(Llines['11']), colorStr+'-.')
    
pylab.xscale('log')    
pylab.yscale('log')
pylab.axis([15, 300, 0.01, 1.0])
pylab.show()    
"""

print 'done'

