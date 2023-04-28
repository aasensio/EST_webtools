from getDDBBdata import getDDBBdata
from photoncount import Photocount
from data_to_csv import toCSV
from data_to_csv import toTransmissions
import fts
import sys
import re

### CONSTANTS
# units
NM_TO_M = 1e-9
RAD_TO_ARCSEC = 206265.

class DatabasePhotoncount:

    def __init__(self, data):
        self.photoncount = Photocount()
        self.data = data
        self.transmissions = []
        self.photoncount.D = 4.2 # telescope diameter in meters
    
    def get_arm_by_wavelength(self, wavelength):
        if (wavelength <= 430):
            return "B"
        elif (wavelength > 430 and wavelength <= 656):
            return "V"
        elif (wavelength > 656 and wavelength <= 854):
            return "R"
        return 'IR'

    def compute_pxscale(self, wavelength):
        pxscale = wavelength * NM_TO_M
        pxscale /= self.photoncount.D
        pxscale *= RAD_TO_ARCSEC
        pxscale /= 2
        return pxscale
    
    def getDL_CLD(self, wavelength):
        arm = self.get_arm_by_wavelength(wavelength)
        if (arm == "B"):
            arm_value = 380
        elif (arm == "V"):
            arm_value = 500
        elif (arm == "R"):
            arm_value = 680
        else:
            arm_value = 1000

        DL = (arm_value * NM_TO_M) / self.photoncount.D * RAD_TO_ARCSEC / 2
        return DL
    
    def compute_pxfactor(self, wavelength, spacres):
        numerator = self.getDL_CLD(wavelength)
        denominator = self.compute_pxscale(wavelength)
        if (spacres != -1): #not dl
            numerator *= spacres
            denominator *= self.getDL_CLD(wavelength)
        pxfactor = numerator / denominator
        return pxfactor
    
    def compute_n_photons(self, nflux, pxfactor, polarimetry):
        n_photons = nflux / 4
        n_photons *= pxfactor ** 2
        if (polarimetry == 1):
            n_photons *= 0.577
        return n_photons
    
    def adjust_nflux(self, nflux, wavelength, specres, int_time):
        nflux *= (wavelength/specres/2 * NM_TO_M) * int_time
        return nflux
    
    def compute_transmissions(self):
        for data in self.data: # if there is no bandpass, snr or integration time. Where wavelength > 1564 errors appear
            #TODO: fixear wavelength si el archivo no se abre
            if (data[2] == "FBI" or data[5] is None or data[6] is None):
                self.transmissions.append([None, 'N/A', None, 'N/A'])
            else:
                #FIXME: asigned bandpass = 1 if there is no bandpass
                if (data[4] is None):
                    data[4] = 1
                lmin = data[0] - (data[4] / 2)
                lmax = data[0] + (data[4] / 2)
                self.photoncount.R = data[3]
                nflux = self.photoncount.compute_nflux(lmin, lmax)
                min_nFlux = nflux.min()
                max_nFlux = nflux.max()
                min_nFlux = self.adjust_nflux(min_nFlux, data[0], data[3], data[6])
                max_nFlux = self.adjust_nflux(max_nFlux, data[0], data[3], data[6])
                
                pxfactor = self.compute_pxfactor(data[0], data[1])
                nphotons_continuum = self.compute_n_photons(max_nFlux, pxfactor, data[7])
                nphotons_core = self.compute_n_photons(min_nFlux, pxfactor, data[7])
                transmission_continuum = data[5] ** 2 / nphotons_continuum
                transmission_core = data[5] ** 2 / nphotons_core
                self.transmissions.append([transmission_continuum, '', transmission_core, ''])
    
    def compute(self):
        self.compute_transmissions()
        toTransmissions(self.transmissions)
        toCSV(self.transmissions)
    
if __name__ == "__main__":
    if (len(sys.argv) < 2 or len(sys.argv) > 3):
        print("Error: wrong number of arguments")
        print("Usage: python3 photoncount.py <path/to/import/file>.csv <path/to/export/file>.csv")
        sys.exit(1)
    match = re.search("\.csv$", sys.argv[1])
    if match is None:
        print("Error: wrong extension on import file")
        print("Usage: python3 photoncount.py <path/to/import/file>.csv")
        sys.exit(1)
    if (len(sys.argv) == 3):
        match = re.search("\.csv$", sys.argv[2])
        if match is None:
            print("Error: wrong extension on export file")
            print("Usage: python3 photoncount.py <path/to/import/file>.csv <path/to/export/file>.csv")
            sys.exit(1)
    data = getDDBBdata(sys.argv[1])
    out = DatabasePhotoncount(data)
    out.compute()