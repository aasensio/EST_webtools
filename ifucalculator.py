import numpy as np
from astropy.io import fits
from scipy.io import readsav
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
import os

"""
Constants and other helpful stuff
"""

IMAGES = {
    "QR": {"name": "Cont500rep.fits", "nx": 1152 * 2, "ny": 1152 * 2, "fact": 10. / 725},  # 1 arcsec ==> 725 km
    "Chromis": {"name": "wb_chromis_3950_05Sep2016_t=54.idlsave", "nx": 1709, "ny": 1122, "fact": 0.0375},
    # arcsec/px = 0.0375 @3950 A
    "Halpha": {"name": "hacore+hawing_24Jun2014_CH_t=312.idlsave", "nx": 833, "ny": 821, "fact": 0.057},  # CRISP Halpha
}
D = 4.2  # M1 diameter in m
RAD_TO_ARCS = 206265
cur_dir = os.path.dirname(__file__)
cur_path = os.path.join(cur_dir, "Images")


class IFUCalculator():

    def __init__(self):
        """
        Computes and displays the FoV for an IFU system depending on the required Spectral and Spatial Resolutions
        re-written in Python by C. Ruiz de Galarreta based on the IDL code from C. Quintero and M. Collados
        """
        """ 
        FoV parameters : 
            Spatial Resolution computed as lambda/D in "*2pix
            FoV
        """
        self.spa_res = None
        self.fov = None

        """ 
        Spectral parameters:
            Spectral Resolution computed as A/2*pix
            Spectral Range that wants to be observed in Angstroms
        """
        self.spec_res = None
        self.spec_range = None

        """ 
        Detector parameters:
            Number of pixels available on each dimension of the detector
            Dual beam configuration: the light is sent to a polarization beamsplitter to produce two orthogonal 
            polarization states. 
            This technique is used for correcting spurious polarization and increase the polarization accuracy.
        """
        self.px_x = None
        self.px_y = None
        self.px_spare = 10  # add an extra number of pixels to be able to differentiate each spectral domain from its
        # adjacent.default is 10 pixels, representing 5 pixels in each side of the spectral domain
        self.beam = None  # default value is single beam so 1

    def computeFPA(self, wavelength, fov, dual_beam=1, simulation=None, px_spare=10):
        # Dual beam
        self.beam = dual_beam

        # fov:
        self.fov = fov

        # Number of pixels required in the spectral dimension for 1 point of the FoV, assuming pixel
        # size in wavelength is half of SpecR
        self.px_spare = px_spare
        spectral_domain = self.spec_range / (self.spec_res / 2) + self.px_spare

        # Size of a squared spatial resolution element in arcsecs^2
        self.IFUparameters(wavelength=wavelength)
        n_2 = (self.spa_res / 2) ** 2

        # Number of required pixels in the Spatial dimension (here Y by convention)
        self.px_y = np.int(self.fov / np.sqrt(n_2))

    def computeFOV(self, wavelength, dual_beam=1, simulation=None, px_spare=10):
        """
        Examples:
            1.
            self.spa_res=0.026        ;[arcsec*2pix]
            self.spec_res=40e-3      ;[A/2*pix]
            self.spec_range=2.      ;[A]
            self.px_x=4000.   ;[pixels]
            self.px_y=4000.   ;[pixels]
        > The message that will appear in your screen will be:
        ifucalculator,0.026,40e-3,2,4000,4000
        The IFU FOV for your configuration is = 4.933 x 4.933 arcsec^2

            2. A second example changing the spatial resolution to 0.1 arcsec
        > ifucalculator,0.100,40e-3,2,4000,4000
        The IFU FOV for your configuration is = 18.97 x 18.97 arcsec^2

            3. third example with a dual beam configuration
        > ifucalculator,0.100,40e-3,2,4000,4000,/Dual_beam
        The IFU FOV for your configuration is = 13.42 x 13.42 arcsec^2

        4. fourth example plotting a context image
        > ifucalculator, 0.056, 60e-3, 5, 8000, 8000,/vis_Ha
        The IFU FOV for your configuration is = 16.80 x 16.80 arcsec^2

         /vis_Ha by /vis_QS or /vis_AR can be exchanged to see different backgrounds
        """
        # Dual beam
        self.beam = dual_beam

        # Number of pixels required in the spectral dimension for 1 point of the FoV, assuming pixel
        # size in wavelength is half of SpecR
        self.px_spare = px_spare
        spectral_domain = self.spec_range / (self.spec_res / 2) + self.px_spare

        # Number of field spectra that  can be packed on the X dimension of the detector
        n_x = np.int(self.px_x / spectral_domain)

        # Total number of field point spectra that can be packed in the Y*X dimensions of the detector
        n_total = self.px_y * n_x / np.int(self.beam)

        # Size of a squared spatial resolution element in arcsecs^2
        self.IFUparameters(wavelength=wavelength)
        n_2 = (self.spa_res / 2) ** 2

        #  The resulting squared FoV * FoV in arcsecs is then:
        self.fov = np.sqrt(n_2 * n_total)

        if simulation == None:
            pass
        else:
            self.plot_simulation()

        return self.fov

    def IFUparameters(self, wavelength):
        # Computes the spatial resolution in " for the given wavelength in nm
        self.spa_res = wavelength * 1e-9 * RAD_TO_ARCS / D

    def load_fits_file(self, file="Cont500rep.fits"):
        fits_path = os.path.join(cur_path, file)
        fits_image = fits.getdata(fits_path, ext=0)
        return fits_image

    def load_idl_save(self, file="hacore+hawing_24Jun2014_CH_t=312.idlsave"):
        idl_path = os.path.join(cur_path, file)
        idl_im = readsav(idl_path)
        # TODO: para Halpha hay dos imagenes del core y de las alas que habrá que sacar
        # idl_im = idl_im["imc"]
        # print(idl_im.shape)
        # plt.figure()
        # plt.imshow(idl_im)
        # plt.show()
        return idl_im

    def plot_simulation(self):
        nx, ny = 1152 * 2, 1152 * 2  # The original simulation has 1152*1152 pixels. They are replicated twice for the plot
        fact = 10. / 725  # The original pixel size is 10 km, here 1" --> 725 kmx
        cx, cy = nx / 2. * fact, ny / 2. * fact

        fits_image = self.load_fits_file()
        plt.figure()
        plt.title('MANCHA Continuum at 500 nm')
        plt.xlabel("X [arcsecs]")
        plt.ylabel("Y [arcsecs]")
        im = plt.imshow(fits_image, cmap='gray', extent=[0, nx * fact, 0, ny * fact], vmin=0.40, vmax=1.3)
        cbar = plt.colorbar(im)
        cbar.set_label("[I/I!dc!n]", rotation=90)

        if self.fov == None:
            pass
        else:
            dim = self.fov / 2.  # in arcsecs
            plt.plot([cx - dim, cx + dim], [cy + dim, cy + dim], color="red")
            plt.plot([cx - dim, cx + dim], [cy - dim, cy - dim], color="red")
            plt.plot([cx - dim, cx - dim], [cy - dim, cy + dim], color="red")
            plt.plot([cx + dim, cx + dim], [cy - dim, cy + dim], color="red")
        plt.show()
