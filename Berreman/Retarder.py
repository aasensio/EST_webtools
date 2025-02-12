"""
Created on Thu Jan 23 12:49:25 2025

@author: esoria
@reviewer: cruiz

Class created from esoria code in order to have a modulable code
"""
import numpy as np
import Berreman
import materials as mat

### CONSTANTS ###

ANGLES = [0, 90, 63.75, 153.75, 0, 90]
pol_waves = [393.3, 396.8, 422.7, 453, 455.4, 460.7, 492.3, 514, 517.3,
                     525, 587.6, 589, 617.3, 630, 656.3, 709, 815.1, 849.8, 854,
                     1083.0, 1565.0]

class Retarder(object):
    def __init__(self):
        self.berreman = Berreman.Berreman()

        #dummy wavelength init:
        self.wavelength = pol_waves[0]

        #quartz crystal orientations:
        self.psi, self.eta, self.xi = 0, 0, 0
        self.eta_angles = self.compute_eta_angles()

        """ optical and mechanical coeffs crom Crystan handbook"""

        #CTEs quartz
        self.cte_o, self.cte_e = 0,0#7.1*1e-6, 13.2*1e-6
        #TOC quartz
        self.dndt_o, self.dndt_e = 13.7*1e-6, 8.9*1e-6#6.5*1e-6, -5.5*1e-6

        #CTEs MgF2
        self.cte_co, self.cte_ce = 0,0 #13.7*1e-6, 8.9*1e-6
        #TOC MgF2
        self.dndt_co, self.dndt_ce = 2.3*1e-6, 1.7*1e-6

        #CTE and TOC infrasil
        self.cte_w, self.dndt_w = 0.55*1e-6, 12.9*1e-6

        """refraction index"""
        ## Mgf2 dummy init
        self.no_c, self.ne_c = None, None

        self.n_oil = 1 #Meadowlark index-matching oil
        #Dummy infrasil init
        self.n_window = None
        ## QUARTZ
        self.no, self.ne = None, None

        self.get_refractiveindex(self.wavelength)

        """thicknessess"""
        ###THICKNESS in nm
        self.PLATE_THICKNESS = [2124900.0, 2099100.0, 2132400.0, 2099100.0, 2124900.0, 2099100.0]
        self.WINDOW_THICKNESS = 10000000.0
        self.OIL_THICNKESS = 10000.0
        self.COATING_THICKNESS = 0*94.4
        ##Optical thicnkess
        self.thkW = None  #window
        self.thkO = None #oil
        self.thkC = None #coating
        self.plate_thickness = None

        self.get_opticalthickness(self.wavelength)

        """field matrix in/out air"""
        self.Fc = self.berreman.fmat([1, 1, 1, 0, 0, 0], 0)  # field matrix for air at normal incidence, n=1, beta=0
        self.Fs = self.berreman.fmat([1, 1, 1, 0, 0, 0], 0)  # field matrix for air output interface

    def compute_eta_angles(self):
        """

        :return: clocking angle for each retarder plate
        """
        return [self.eta + np.radians(angle) for angle in ANGLES]

    def get_refractiveindex(self, wl, DT=0.):
        """

        :param wl: wavelength in nanometers
        :return:
        """
        inx2 = mat.mgf2(wl)
        self.no_c, self.ne_c = inx2[0, 0], inx2[0, 1]
        self.no_c *= (1+self.dndt_co*DT)
        self.ne_c *= (1+self.dndt_ce*DT)

        self.n_window = mat.fusedsilica(wl)
        self.n_window *= (1+self.dndt_w*DT)

        self.no, self.ne = mat.quartz_retarder(wl)
        self.no *= (1+self.dndt_o*DT)
        self.ne *= (1+self.dndt_e*DT)

    def get_opticalthickness(self, wavelength, DT=0):
        """
        Computes he opti
        :param wavelength:
        :return: updates the variables for optical thickness
        """
        self.thkW = self.compute_optical_thickness(self.WINDOW_THICKNESS, wavelength, CTE=self.cte_w, DT=DT) #window
        self.thkO = self.compute_optical_thickness(self.OIL_THICNKESS, wavelength, DT=DT) #oil
        self.thkC = self.compute_optical_thickness(self.COATING_THICKNESS, wavelength, CTE=self.cte_co, DT=DT) #coating
        self.plate_thickness = self.compute_optical_thickness(self.PLATE_THICKNESS, wavelength, CTE=self.cte_o,DT=DT)

    def compute_optical_thickness(self, thickness, wavelength, CTE=0, DT=0):
        """

        :param thickness: mechanical thicnkess
        :param wavelength: wavelencth in nm
        :return: optical thickness for the given wavelength
        """
        if isinstance(thickness, list):
            return [(t*(1+CTE*DT))/ wavelength for t in thickness]
        else:
            return thickness*(1+CTE*DT)/wavelength

    def create_layer_stack(self, bonding=None, window=True):
        """
        Creates a zero order plate made of a 6 layer of stacked quartz
        :param bonding: the bonding method can be oil or optical contact
        :param window: f an input/output window is applied, works only for oil bonding
        :return: the stack layer
        """
        layers = []

        if bonding == "oil":
            if window: # case with infrasil in/out window
                layers += [
                    [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC],  # Input Coating
                    [self.n_window, self.n_window, self.n_window, 0, 0, 0, self.thkW],  # Input Window
                    [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC],  # Coating after Window
                    [self.n_oil, self.n_oil, self.n_oil, 0, 0, 0, self.thkO]  # Oil Layer
                ]
            # Main stack: (Coating → Quartz → Coating → Oil) x 5
            for thk, eta in zip(self.plate_thickness[:-1], self.eta_angles[:-1]):
                layers += [
                    [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC],  # Coating
                    [self.no, self.ne, self.no, eta, self.psi, self.xi, thk],  # Quartz
                    [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC],  # Coating
                    [self.n_oil, self.n_oil, self.n_oil, 0, 0, 0, self.thkO]  # Oil
                ]
            # Final coating → quartz → coating
            layers += [
                [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC],  # Coating
                [self.no, self.ne, self.no, self.eta_angles[-1], self.psi, self.xi, self.plate_thickness[-1]],
                # Final Quartz
                [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC]  # Coating
            ]
            # Add output window if included
            if window:
                layers += [
                    [self.n_oil, self.n_oil, self.n_oil, 0, 0, 0, self.thkO],  # Oil
                    [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC],  # Coating
                    [self.n_window, self.n_window, self.n_window, 0, 0, 0, self.thkW],  # Output Window
                    [self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC],  # Output Coating
                ]
        elif bonding == "contact": #No window
            layers += [[self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC]] # top coating
            #Uncoated quartz stack
            for thk, eta in zip(self.plate_thickness, self.eta_angles):
                layers += [[self.no, self.ne, self.no, eta, self.psi, self.xi, thk]]
            layers += [[self.ne_c, self.no_c, self.no_c, 0, 0, 0, self.thkC]]  # bottom coating

        return np.asarray(layers)

    def compute_mueller_stack(self, wavelength, bonding=None, window=True, DT=0):
        """

        :param wavelength: wavelength in nm
        :param bonding: the bonding method can be either optical contact or oil
        :param window: only needed for oil bonding
        :return: Transmitted mueller matrix for the stack at a given wavelength
        """
        #update refraction index values
        self.get_refractiveindex(wavelength, DT=DT)

        #updated optical thicness
        self.get_opticalthickness(wavelength, DT=DT)

        #create layer stack
        layer = self.create_layer_stack(bonding=bonding, window=window)
        #compute mueller matrix elements
        M = self.berreman.cmat(layer, 0.0) # characteristics matrix
        MMr, MMt, Sr, St = self.berreman.cmmat(self.Fc, M, self.Fs)

        #Normalize the transmission matrix
        MMT_00 = MMt[0, 0]
        MMt /= MMT_00
        MMt[0, 0] = MMT_00
        #FIXME:throws an inf error when computing with an oil+window config
        return MMt

    def montecarlo(self, s_mean, s_std, N=1000):
        samples = np.random.normal(s_mean, s_std, N)
        return samples




