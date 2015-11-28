from numpy import zeros, fabs, arange, array_equal, exp, ones
from numpy import linalg, eye, dot
from IPython.core.debugger import Tracer

import read_ei
import read_cr
import read_levels
from scipy.constants import c,h,Boltzmann
from Leiden_ISM.ismUtils import planckOccupation as ng

import temperature

kb = Boltzmann

def reduce_vj_repr(en, a_eins, cr, T,  ini, fin, vj_unique,
                   debug=False):
    """.. todo:: add doc"""

    # number of v and j levels respectively
    nv, nj =  en.shape

    # linear indices for all the energy levels
    lind = arange(en.size)

    # the linear indices indexed as a 2D array
    # lind2d[v,j] = the linear index
    lind2d = lind.copy().reshape(en.shape)

    # convert the energies from a v,j representation
    # to a linear representation
    en_l = en.flatten()
    if debug is True:
        assert en[:,:] - en_l[lind2d[:,:]]

    # convert the representation of the einstein
    # coefficients from a 4D representation v,j,v',j'
    # to a 2D representation i,j
    a_eins_l = a_eins.reshape(en.size, en.size)
    # added for checking the reshaping
    if debug is True:
        for v in arange(nv):
            for j in arange(nj):
                for vp in arange(nv):
                    for jp in arange(nj):
                        a_l = a_eins_l[lind2d[v,j], lind2d[vp,jp]]
                        a = a_eins[v,j,vp,jp]
                        assert a_l == a

    # initializing the linearized version of the
    # the collisional reaction rates and 
    # initial and final quantum numbers; the  
    # reshape is done according to the total 
    # number of rovibrational levels of H2 
    cr_l = zeros((en.size, en.size, T.size), 'f8')
    ini_l = zeros(ini.shape[1], 'i')
    fin_l = zeros(fin.shape[1], 'i')
    
    for i, ((v,j), (vp,jp)) in enumerate(zip(ini.T, fin.T)):
        cr_l[lind2d[v,j], lind2d[vp,jp], :] = cr[v, j, vp, jp,:]
        ini_l[i] = lind2d[v,j]
        fin_l[i] = lind2d[vp,jp]
    vj_unique_l = lind2d[vj_unique[0], vj_unique[1]]
    # added for checking the reshaping
    # .. todo:: think of a better test
    # .. todo:: see how ini_l and fin_l can be checked
    if debug is True:
        for (v,j), (vp,jp) in zip(ini.T, fin.T):
            T1 = cr_l[lind2d[v,j], lind2d[vp,jp], :]
            T2 = cr[v, j, vp, jp, :]
            assert array_equal(T1, T2)
    

    # computing the degeneracies
    g = vj_unique[1]*2 + 1

    return en_l, a_eins_l, cr_l, ini_l, fin_l, vj_unique_l, g


#def computeRateMatrix(pNH3, Tkin, nc):
def computeRateMatrix(en, a_eins, cr, ini, fin, unique_col, g, tkin, nc):
    """compute the matrix of transition rates"""

    ###########################################################
    # constructing the matrix
    ###########################################################
    def fill_K_matrix(cr=None,
                      levels=None,
                      collider_density=None,
                      unique_col=None,
                      g=None,
                      ini=None,
                      fin=None, 
                      tkin=None):
        """fill the kij matrix from the collsion rates
        cr:  matrix containing the coolisional
             rates in the linearized form (58,58,50)
        levels: energies for all the rovibrational 
                state of H2 (480 from UCLA database)
        collider_density: density of the coolider for which
                          the collisional reaction rates are
                          given (H in Francois's case)
        unique_col: couples (v,j) for which the collisional  
                    reaction rates are provided
                    (58 in Francois's case)
        g: degeneracies of the levels for which the collisional 
           reaction rates are provided
        ini: couple of all possible initial (v,j) in the provided
             table of reaction rates
        fin: couple of all possible final (v',j') in the provided
             table of reaction rates
        tkin: array containing the temperature at which the reaction
              rates are provided
        .. todo:: improve documentation
        """

        # 1) get the rate coefficient for level i -> i' at
        # temperature T. The collision rate matrix should
        # look like this
        #
        #    k = rate_coeff_for_transition_cr[i, i', T]
        # 
        # 2) get the energy difference between i -> i'
        #
        #    e = fabs(E[i] - E[i']
        #
        # 3) fill the k matrix by repeating 1 for the direct transitions
        #    i.e. i -> i'. and use the detailed balance to
        #    get the rate coeff for the reverse transtion i' -> i

        n = unique_col.size
        K = zeros((n, n), 'f8')

        # mapping the level number to the entry column/row
        # in the K matrix
        ired = zeros(unique_col.max() + 1, 'i4')
        for i, iu in enumerate(unique_col):
            ired[iu] = i

        # filling the K matrix
        for i, ip in zip(ini, fin):

            # indicies of the levels in the matrix K
            ir, irp = ired[i], ired[ip]

            # the difference between the two energy levels
            dE = fabs(levels[i] - levels[ip])

            # the collision rate at tkin[0]. we are useing ir and ip in
            # indexing g since g has the same length as the unqiue levels
            tkin_ind = 40
            criip = cr[i,ip, tkin_ind]
            criip_rev = cr[ip,i, tkin_ind] * exp( -dE / (kb*tkin[tkin_ind]) ) * g[ir]/g[irp]
            
            K[ir,irp] = criip
            K[irp,ir] = criip_rev
        
        return K


    def fill_AP_matrix(a_eins=None,
                       levels=None,
                       unique_col=None,
                       ini=None,
                       fin=None,
                       tkin=None):
        """fill the (A prime)_ij matrix from the lambda radiative transitions
        """
        tcmb = tkin
        n = unique_col.size
        AP = zeros( (n, n), dtype = 'f8')

        # mapping the level number to the entry column/row
        # in the K matrix
        ired = zeros(unique_col.max() + 1, 'i4')
        for i, iu in enumerate(unique_col):
            ired[iu] = i

        # filling the AP matrix
        for i, ip in zip(ini, fin):

            # indicies of the levels in the matrix K
            ir, irp = ired[i], ired[ip]

            # the difference between the two energy levels
            dE = fabs(levels[i] - levels[ip])
            
            nu = dE / h # freq in Hz
 
            # the einstein coefficients. we are useing ir and ip in
            # indexing g since g has the same length as the unique levels   
            tcmb_index = 40
            APiip = (1.0 + ng(h, nu, kb, tcmb[tcmb_index]))*a_eins[i,ip]
             
            AP[ir,irp] = APiip
        return AP
  
    def fill_ABS_matrix(a_eins=None,
                       levels=None,
                       unique_col=None,
                       g = None,
                       ini=None,
                       fin=None,
                       tkin=None):
        """fill the Aij matrix for absorbtion transitions from the lambda
        radiative transitions"""
        tcmb = tkin
        n = unique_col.size
        ABS = zeros( (n, n), dtype = 'f8')

        # mapping the level number to the entry column/row
        # in the K matrix
        ired = zeros(unique_col.max() + 1, 'i4')
        for i, iu in enumerate(unique_col):
            ired[iu] = i

        # filling the AP matrix
        for i, ip in zip(ini, fin):

            # indicies of the levels in the matrix K
            ir, irp = ired[i], ired[ip]

            # the difference between the two energy levels
            dE = fabs(levels[i] - levels[ip])
            nu = dE / h # freq in Hz
 
            # the einstein coefficients. we are useing ir and ip in
            # indexing g since g has the same length as the unique levels   
            tcmb_index = 40
            ABSiip = (g[ir]/g[irp])*ng(h, nu, kb, tcmb[tcmb_index])*a_eins[i,ip] 
            ABS[ir,irp] = ABSiip
        
        return ABS
    
    def fill_E_matrix(unique_col = None):
        """This is a matrix full of ones (it is a utility matrix"""
        n = unique_col.size
        return ones((n, n))


    n = unique_col.size

    k_mat = fill_K_matrix(cr=cr,
                          levels=en,
                          collider_density=nc,
                          unique_col=unique_col,
                          g=g,
                          ini=ini,
                          fin=fin,
                          tkin=tkin)

    ap_mat = fill_AP_matrix(a_eins=a_eins,
                            levels=en,
                            unique_col=unique_col,
                            ini=ini,
                            fin=fin,
                            tkin=tkin)

    abs_mat = fill_ABS_matrix(a_eins=a_eins,
                              levels=en,
                              unique_col=unique_col,
                              g = g,
                              ini=ini,
                              fin=fin,
                              tkin=tkin)

    E = fill_E_matrix(unique_col = unique_col)

    F = nc * k_mat + ap_mat + abs_mat.T
    diag = eye(n)*dot(F, E)[:, 1]
    offdiag = F.T
    full = -diag + offdiag
    
    return full


def solveEquilibrium(full):
    """solve for the equlibrium population densities.
    
    i.e solving A.x = b
    where full = A
    """

    n = full.shape[0]
    
    # solving directly
    # replacing the first row with the conservation equation
    dndt = zeros((n, 1), 'f8')
    full[0,:] = 1.0
    dndt[0]   = 1.0

    # solving the system A.x = b
    # before solving, we will devide each row by the diagonal
    A, b = full, dndt
    for i in arange(n):
        A[i,:] = A[i,:]/A[i,i]
        x = linalg.solve(A, b)

    #print x.T

    # the fractional population density
    f = x
    return f

#def Tg2Tr(Tg):
#    '''Given the array of kinetic temperatures, it returns the
#     array with the corresponding radiation temperatures (through
#     the redshift parameter z)
#     '''
#     return Tr



def coolingFunction(x, en, eins, T, ini, fin, unique_col):
    cf = zeros(1,'f8')
    ired = zeros(unique_col.max() + 1, 'i4')
    for i, iu in enumerate(unique_col):
        ired[iu] = i
        tkin_index = 40

    for i, ip in zip(ini, fin):
        # indicies of the levels in the matrix K
            ir, irp = ired[i], ired[ip]

            # the difference between the two energy levels:
            dE = fabs(en[i] - en[ip])

            # the einstein coefficient connecting i to ip:
            Ai_ip = eins[i,ip]

            #the population of the upper level:
            chi = x[i]
            cf =+  Ai_ip*dE*x[i]
    return cf