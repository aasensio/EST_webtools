import os.path
from scipy.io import readsav

DATAFILES = [
    './fts_disk_center.idlsave',
'./ftscnt_full.sav']

CLIGHT = 2.99792458e8  # speed of light [m/s]
AA_TO_M = 1e-10
CM_TO_M = 1e-2
ERG_TO_WATT = 1e-7

class fts(object):

    def __init__(self):
        self.ll = None
        self.ii = None
        self.cc = None
        self.nu = None

        self.datafile = './fts_disk_center.idlsave'

    def get_atlas(self, atlas=0):

        # The atlas parameter defines if the reference atlas is the original from J. Leenarts, expressed in watt/(cm2 ster AA)
        # as emitted at solar surface
        # Or the new atlas exported by M. Collados allowing for going beyong in the IR spectrum, and expressed
        # in erg/s/(cm2sterAA)
        self.datafile = DATAFILES[atlas]

        try:
            t = readsav(self.datafile)
        except FileNotFoundError:
            print("Error Atlas not found", DATAFILES[1])
            atlas = 0
            self.datafile = DATAFILES[0]
        finally:
            print("Opening Atlas", DATAFILES[atlas])

        t = readsav(self.datafile)

        if atlas == 0:
            # convert to J s-1 m-2 m-1 sr-1
            self.ll = t['ftswav'] * AA_TO_M
            self.nu = CLIGHT / self.ll
            self.ii = t['ftsint'] * CM_TO_M ** (-2) * AA_TO_M ** (-1)  # from from W /( cm2 ster AA) to  W / (m2 ster m)
            self.cc = t['ftscnt'] * CM_TO_M ** (-2) * AA_TO_M ** (-1)

        else:
            ### for ftscnt_full.sav
            ## lam: wavelength in A  with 2 mA steps. IR data has been interpolated as the step was 4 mA.
            ## ftsnorm:spectrum normalized to the continuum
            ## cont: value of the continuum at all wavelengths.
            # The spectro in erg/s/cm2/A/sr) is the product of ftsnorm and cont.
            self.ll = t["lam"] * AA_TO_M
            cont = t["cont"]
            ftsnorm = t["ftsnorm"]
            self.ii = cont * ftsnorm * ERG_TO_WATT * CM_TO_M ** (-2) * AA_TO_M ** (-1)





       
