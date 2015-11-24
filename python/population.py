import read_ei
import read_cr
import read_levels


#def computeRateMatrix(pNH3, Tkin, nc):
def computeRateMatrix(nc):
    """compute the matrix of transition rates"""

    n = pNH3.nlevels

    levels = read_levels.read_levels('Read/H2Xvjlevels.cs')
    transColl, T, ini, fin = read_cr.read_coeff("Read/Rates_H_H2.dat")
    transRad = read_ei.read_einstein()

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