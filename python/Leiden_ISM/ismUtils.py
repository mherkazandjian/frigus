# utility function useful for ISM stuff
import numpy as np
from scipy import integrate
import numpy

def Av2NH(Av, Z):
    """Computes the column density of hydrogen nuclei given an Av using the formula
    in paper1 (add ref)
    """
    return (Av*1.87e21)/Z
def NH2Av(NH, Z):
    """The inverse of :data:`Av2NH`"""
    return NH*(1.0/1.87e21)*Z

    
def AvToLength(Av, nGas, Z):
    """ convert the input visual extinction value to the length in cm (L) of the the NH 
        column (not the actual column density) this is the inverse of Eq-4 
        in paper01 where  :math:`N_H = L.n_{gas}` 
        
        :param Av: The Av (in mag units).
        :type Av: numpy.float64
        :param nGas: The gas density in :math:`cm^{-3}`.
        :type nGas: numpy.float64
        :param Z: The metllicity in terms of solar metallicity (:math:`Z_{\odot}`).
        :type Z: numpy.float64
    """
    return (Av*1.87e21)/(nGas*Z)

def LengthToAv(NH, nGas, Z):
    """ convert the column density to visual extinction value in Av. This is
        Eq-4 in paper01 implemented in fucntion :data:`LengthToAv`
        
    :param NH: The H column density (in :math:`cm^{-2}`).
    :type NH: numpy.float64
    :param nGas: The gas density in :math:`cm^{-3}`.
    :type nGas: numpy.float64
    :param Z: The metllicity in terms of solar metallicity (:math:`Z_{\odot}`).
    :type Z: numpy.float64
    """
    return NH*nGas*Z / 1.87e21

def getSlabThicknessFromAv(slabsAv, nGas, Z):
    """ given the Av of the sub-slabs of a plane parallel model, 
        it computes the thickness of each sub-slab of the descretized
        model. All the parameters are of the same type as :data:`AvToLength`
        except for :
        
        :param slabsAv: An array holding the Av (in mag units) for the staring positions of the sub-slabs.
        :type slabsAv: numpy.ndarry([], dtype = numpy.float64)
        :returns: The thickness of each slab in cm
        :rtype: numpy.ndarry([], dtype = numpy.float64)
     """
    dAvSlabs = slabsAv[1::] - slabsAv[0:-1]
    dxSlabs  = AvToLength( dAvSlabs, nGas, Z)
    return dxSlabs
    
def planckOccupation(h, nu, kb, T):
    """computes the planck function for input parameters.
       :keywords: beta, nu,  
    """
    x = h * nu / (kb * T)
    return 1.0 / (np.exp(x) - 1.0)

def ortho_para_abundance_at_eq(tkin, xH2):
    """compute the ortho to para abundance (relative to xH2) at thermodynamic equilibirum given a kinetic temperature.
    i.e here, they add up to """
    rop = 9.0*np.exp(-170.6/tkin)
    
    xoH2 = xH2 / (1.0 + 1.0 / rop)
    xpH2 = xH2 - xoH2
    
    return xoH2, xpH2


def startburst_gamma_mech_rate(SFR=None, n_PDR=None, d_PDR=None, d_SB=None, E_SN=None, eta=None):
    '''computes the mechanical heating rate in a startburst reigon given :
    
       - star formation rate in Msun/yr (SFR) 
       - number density of PDRs per pc^3
       - the average diameter of each PDR in pc
       - energy of each super nova (E_SN) event = 1e51 egrs
       - SN heating efficiency (the efficiency of turbulent heating due to the SN shock absorbed by the ISM)
    '''

    pc2cm = 3.08e18
    yr2sec = 365.25 * 24.0 * 3600.0
    
    def f_IMF(x):
        return x**(-2.35)
    
    #relation between SFR and SNR, k = SNR / SFR
    mass_SN_rng = numpy.linspace(8.0, 50.0, 1000.0)
    mass_SF_rng = numpy.linspace(0.1, 125.0, 10000.0)
    k = integrate.simps(f_IMF(mass_SN_rng), mass_SN_rng) / integrate.simps(mass_SF_rng*f_IMF(mass_SF_rng), mass_SF_rng) 
    
    #the supernova rate
    SNR = SFR * k
    
    #volume of each PDR (in cm^3)
    V_PDR = (4.0*numpy.pi/3.0)*(0.5*d_PDR*pc2cm)**3.0
    
    #volume of starburst reigon (in pc^3, yes, it is pc^3 not cm^3)
    V_SB = (4.0*numpy.pi/3.0)*(0.5*d_SB)**3.0 
    
    #number of PDRs in the SB
    N_PDR = V_SB * n_PDR
    
    #energy input absorbed by the PDRs due to all the supernovea (erg / sec)
    E_SN_absorbed = SNR * E_SN * eta / yr2sec
    
    #total volume of PDRs in the starburst 
    V_total_PDRs_in_SB = N_PDR * V_PDR 
    
    #mechanical heating rate in erg/cm^3/s
    #mean_gamma_mech_per_unit_volume_in_PDRs = E_SN_absorbed / (total_NR * E_SN * eta /) 
    gamma_mech = E_SN_absorbed / V_total_PDRs_in_SB 
    
    print gamma_mech
    
    
def LN_dispersion_to_mach_num(sigma):
    '''computes and estimate of the mach number given the width of the Log normal denisty
    distribution'''
    
    return numpy.sqrt( (numpy.exp(sigma**2.0) - 1.0)*(4.0/3.0))
