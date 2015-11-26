from numpy import zeros, fabs, arange
from IPython.core.debugger import Tracer

import read_ei
import read_cr
import read_levels
from scipy.constants import c,h,Boltzmann

kB = Boltzmann

def reduce_vj_repr(en, a_eins, cr, ini, fin, vj_unique):

    nv, nj =  en.shape

    # linear indices for all the energy levels
    lind = arange(en.size)

    # the linear indices indexed as a 2D array
    # lind2d[v,j] = the linear index
    lind2d = lind.copy().reshape(en.shape)

    # convert the energies from a v,j representation
    # to a linear reporesentation
    en_l = en.flatten()
    #assert en[:,:] - en_l[lind2d[:,:]]

    # convert the representation of the einstein
    # coefficients from a 4D representation v,j,v',j'
    # to a 2D representation i,j
    a_eins_l = a_eins.reshape(en.size, en.size)
    #for v in arange(nv):
    #    for j in arange(nj):
    #        for vp in arange(nv):
    #            for jp in arange(nj):
    #                a_l = a_eins_l[lind2d[v,j], lind2d[vp,jp]]
    #                a = a_eins[v,j,vp,jp]
    #                assert a_l == a + 0.1

    cr_l = zeros((en.size, en.size, T.size), 'f8')
    for (v,j), (vp,jp) in zip(ini.T, fin.T):
        cr_l[lind2d[v,j], lind2d[vp,jp], :] = cr[v,j,vp,jp,:]
    ### add the check here ###
    # asdasdasd
    ##########################
    Tracer()()
    

    return en_l, cr_l, a_eins_l, ini_l, fin_l, vj_unique_l


#def computeRateMatrix(pNH3, Tkin, nc):
def computeRateMatrix():
    """compute the matrix of transition rates"""
    rovib_levels = read_levels.read_levels("Read/H2Xvjlevels.cs")
    transColl, T, ini, fin, vj_unique = read_cr.read_coeff("Read/Rates_H_H2.dat")

    A = read_ei.read_einstein()

    ###########################################################
    # constructing the matrix
    ###########################################################
    def fill_K_matrix(cr=None,
                      levels=None,
                      collider_density=None, tkin=None):
        """fill the kij matrix from the collsion rates
        
        .. todo:: improve documentation
        """

        # 1) get the rate coefficient for (v,j) -> (v',j')
        # at temperature T
        #
        #    k = rate_coeff_for_transition_cr[v,j, v', j', T]
        # 
        # 2) get the energy difference between (v,j) -> (v',j')
        #
        #    e = fabs(E[v,j] - E[v',j']
        #
        # 3) fill the k matrix by repeating 1 for the direct transitions
        #    i.e. (v,j) -> (v',j'). and use the detailed balance to
        #    get the rate coeff for the reverse transtion (v',j') -> (v,j)

        n = vj_unique.shape[1]
        K = zeros( (n, n), 'f8')

        for (v,j), (vp,jp) in zip(ini.T, fin.T):

            # the difference between the two energy levels
            dE = fabs(levels[v, j] - levels[vp, jp])

            #attaching a label to each unique level
        for i, couple in vj_unique.T:
            couple_label[i] = i

        #for trans in transColl:
            

        """
        for trans in transColl:
            # id1_trans =
            # id2_trans =
            # [vj_unique[0][:],vj_unique[1][:]][v[trans]][j[trans]]
            u = ini[:,trans]
            l = fin[:,trans]
            gu = 2*u[1]+1
            gl = 2*l[1]+1
            v = ini[:,trans][0]
        #cr[ini[:,1][0],ini[:,1][1],fin[:,1][0],fin[:,1][1]]

            # difference in the energy in K
            dE = abs(levels[v,j] - levels[vp,jp])
            K[u,l] = transColl[ini[:,trans][0],ini[:,trans][1],fin[:,trans][0],fin[:,trans][1]]
        #            K[u, l] = trans['rc'](Tkin)
        #            K[l, u] = (float64(gu) / float64(gl)) * K[u,l] * exp(-dE / Tkin)
        """
        return K
        
    k_mat = fill_K_matrix(cr=transColl,
                          levels=rovib_levels,
                          collider_density=2.5,
                          tkin=400.0)
    asdasdasd
    # def fill_AP_matrix():
    #     """fill the (A prime)_ij matrix from the lambda radiative transitions
    #     database"""
    #     AP = zeros( (n, n), dtype = float64)
    #
    #     for trans in transRad:
    #         u  = trans['u']; l = trans['l']
    #         dE = abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
    #
    #         nu = dE*kBoltz / hPlank # freq in Hz
    #
    #         AP[u, l] = (1.0 + ng(hPlank, nu, kBoltz, Tcmb))*trans['A']
    #
    #     return AP
    #
    # def fill_ABS_matrix():
    #     """fill the Aij matrix for absorbtion transitions from the lambda
    #     radiative transitions database"""
    #     ABS = zeros( (n, n), dtype = float64)
    #
    #     for trans in transRad:
    #         u  = trans['u']; l = trans['l']
    #         dE = abs( levels[u]['E'] - levels[l]['E'] ) # energy in K
    #         gu, gl = float64(levels[u]['g']), float64(levels[l]['g'])
    #
    #         nu = dE*kBoltz / hPlank # freq in Hz
    #
    #         ABS[u, l] = (gu/gl)*ng(hPlank, nu, kBoltz, Tcmb)*trans['A']
    #
    #     return ABS
    #
    # def fill_E_matrix():
    #     """This is a matrix full of ones (it is a utility matrix"""
    #     return ones((n,n))

    K = fill_K_matrix()
    AP = fill_AP_matrix()
    ABS = fill_ABS_matrix()
    E = fill_E_matrix()

    F = nc * K + AP + ABS.T
    diag = eye(n)*dot(F, E)[:, 1]
    offdiag = F.T

    full = -diag + offdiag

    return full


# def solveEquilibrium(pNH3, full):
#     """solve for the equlibrium population densities"""
#     n = pNH3.nlevels
#
#     # solving directly
#     #replacing the first row with the conservation equation
#     dndt = zeros((n,1))
#     full[0,:] = 1.0
#     dndt[0]   = 1.0
#
#     A = full
#     b = dndt
#     #solving the system A.x = b
#     #before solving, we will devide each row by the diagonal
#     for i in arange(n):
#         A[i,:] = A[i,:]/A[i,i]
#     x = linalg.solve(A, b)
#
#     #print x.T
#
#     # the fractional population density
#     f = x
#     return f
