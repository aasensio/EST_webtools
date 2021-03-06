import numpy as np
from scipy.signal import fftconvolve
import fts
import matplotlib.pyplot as pl


class Photocount(object):

    def __init__(self):
        """
        Computes and displays the integration to reach a desired level
        of S/N per diffraction limited spatial pixel and given spectral
        resolution.
        
        written by Jorrit Leenaarts
        """

        # units
        self.hplanck = 6.626e-34 # m**2 kg /s
        self.clight = 2.99792458e8 #m /s
        self.eV = 1.602e-19 # J
        self.mAA_to_m = 1e-13
        self.nm_to_m = 1e-9
        self.m_to_nm = 1e9
        self.nm_to_pm = 1e3
        self.rad_to_arcsec = 206265.
        self.km_to_m = 1e3
        self.m_to_arcsec = 1./(7.25e5)

        self.atlas = fts.fts()
        self.al, self.at = self.get_atmostrans()

        self.D = 4.0
        self.lmin = 854.0
        self.lmax = 855.0
        self.polarimetry = 0
        self.R = 8e4
        self.T = 0.1
        self.SN = 1e3
        self.v = 7.0
        self.binning = 1.0
        
        self.resmin = self.spatres(self.lmin * self.nm_to_m)
        self.resmax = self.spatres(self.lmax * self.nm_to_m)

        self.sresmin = self.specres(self.lmin * self.nm_to_m)
        self.sresmax = self.specres(self.lmax * self.nm_to_m)

    def get_atmostrans(self, M = 2.0):
        """
        Reads atmospheric transmission at La Palma from "D. L. King,
        Atmospheric Extinction at the Roque de los Muchachos
        Observatory, La Palma, RGO/La Palma technical note no 31"
        available at http://www.ing.iac.es/Astronomy/observing/manuals/man_tn.html
    
        M sets the elevation. M=2 means 30 degrees, zenith is M=1
        """

        a = np.loadtxt('atmostrans.dat')
        ll = a[:,0] *1e-10 # from AA to m
        T = a[:,1]
        TT = 10.**(-0.4 * T * M )
        return ll, TT

    def get_filter(self, ll, fwhm):
        #
        # just a normalized gaussian
        #
        llc= ll[int(len(ll)/2)]
        sigma = fwhm / (2.0 * (2.0 * np.log(2.0))**0.5)
        trans = np.exp(-0.5 * ((ll - llc) / sigma)**2)
        trans /= trans.sum()
        return trans

    def spatres(self, l):
        out = 1.22 * l / self.D * self.rad_to_arcsec / 2.
        return out

    def specres(self, l):
        out = l / self.R / 2.0 * self.nm_to_pm
        return out

    def compute(self):
        if (self.polarimetry == 0):
            pfac = 1.0
        else:
            pfac = 1.0 / 3.0

        if (self.lmin >= self.lmax):
            self.lmin = self.lmax - 1.0

        self.ran = (self.lmin, self.lmax)
        lmin = self.lmin - 0.1
        lmax = self.lmax + 0.1

        # cut out selected  wavelength range
        iw = np.where( (self.atlas.ll > lmin * self.nm_to_m) &
                      (self.atlas.ll < lmax * self.nm_to_m))[0]
        self.Ilambda = self.atlas.ii[iw]
        self.ll = self.atlas.ll[iw]

        # interpolate atmospheric transmission
        atrans = np.interp(self.ll, self.al, self.at) 

        # wavelength spacing
        dl = self.ll / self.R

        # spectrally convolved spectrum, assume wavelength spacing is
        # equidistant.  the spacing changes slowly over the entire
        # wavelength range so this is approx OK.
        iw = int( len(self.ll) / 2)
        dl2 = self.ll[iw+1] - self.ll[iw]
        ff = self.get_filter(self.ll, self.ll[iw] / self.R)
        self.Ilambdac = fftconvolve(self.Ilambda, ff, mode = 'same')

        # energy per photon [J]
        ephot = self.hplanck * self.clight / self.ll

        # nphot = Ilambda / ephot * A              * T * (delta omega)**2               * (delta lambda)
        #       = Ilambda / ephot * f0 * pi * R**2 * T * 1.22**2 * lambda**2 / (4 R**2) * (delta lambda)
        #       = Ilambda / ephot * f0 * pi        * T * 1.22**2 * lambda**2 / 4        * (delta lambda)

        # f0 is the fraction of the aperture that is unobscured. I
        # assume 1.0 here

        # number of photons per second per m per diffraction limited
        # spatial resolution element of the light entering the
        # telescope
        nflux = self.Ilambda / ephot * np.pi * atrans * 1.22**2 / 4.0 * self.ll**2

        # number of photons per second per spatial pixel
        nflux /=  4.0 
         # per spectral pixel
        nflux *= (dl / 2.0)
         # factor in spatial binning
        nflux *= self.binning**2
         # and factor in polarimetry
        nflux *= pfac
 
        nfluxideal = np.copy(nflux) # ideal telescope
        nflux *= self.T # for given total transmission

        # time to reach desired S/N, computed from: S/N = sqrt(nflux * t)
        self.t = self.SN**2 / nflux
        self.tideal = self.SN**2 / nfluxideal

        # display spatial resolution
        self.resmin = self.spatres(self.binning * lmin)
        self.resmax = self.spatres(self.binning * lmax)

        self.sresmin = self.specres(lmin)
        self.sresmax = self.specres(lmax)

        # compute Alex Feller ideal dt,dx
        phi = self.Ilambda / ephot * np.pi / 4.0 * self.D**2 * atrans * self.T * pfac # photons/ (s ster m)
        phi = phi * self.nm_to_m / self.rad_to_arcsec**2 # photons/ (s arcsec**2 nm)
        self.dt = self.SN**2 / phi / (dl/self.nm_to_m/2.0) / (self.v* self.km_to_m*self.m_to_arcsec)**2
        self.dt = self.dt**(1./3.)
        self.dx = (self.v* self.km_to_m*self.m_to_arcsec) * self.dt 

    def plot(self):
        f, ax = pl.subplots(ncols=2, nrows=2, figsize=(10,10))

        xax = self.ll * self.m_to_nm

        ax[0,0].plot(xax, self.Ilambda, label = 'atlas')
        ax[0,0].plot(xax, self.Ilambdac, label = 'smeared to R')
        ax[0,0].legend(loc='best')
        ax[0,0].grid(True)
        ax[0,0].set_xticklabels([])
        ax[0,0].set(xlim=self.ran, ylabel=r'$I_\lambda$ [W m$^{-2}$ m$^{-1}$ sr$^{-1}$]', title='Spectrum')
        
        ax[0,1].plot(xax, self.dx)
        ax[0,1].grid(True)
        ax[0,1].set_xticklabels([])
        ax[0,0].set(xlim=self.ran, ylabel=r'$\Delta x$ [arcsec]', title='optimal pixel size for given signal speed')

        ax[1,0].plot(xax, self.t, label='with given transmission')
        ax[1,0].plot(xax, self.tideal, label='perfect telescope')
        ax[1,0].grid(True)
        ax[1,0].legend(loc='best')
        ax[1,0].set(xlim=self.ran, xlabel=r'$\lambda$ (nm)', ylabel=r'$\Delta t$ [s]', title='integration time for given spatial and spectral pixel size')

        ax[1,1].plot(xax, self.dt)
        ax[1,1].grid(True)
        ax[1,1].set(xlim=self.ran, xlabel=r'$\lambda$ (nm)', ylabel=r'$\Delta t$ [s]', title='optimal integration time for given signal speed')

        pl.show()

    def set_properties(self, d):

        self.D = d['D']
        self.lmin = d['lmin']
        self.lmax = d['lmax']
        self.polarimetry = d['polarimetry']
        self.R = d['R']
        self.T = d['T']
        self.SN = d['SN']
        self.v = d['v']
        self.binning = d['binning']

if (__name__ == '__main__'):
    out = Photocount()
    out.compute()
    out.plot()