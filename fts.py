import numpy as np
from scipy.io import readsav

class fts:

    ll = None
    ii = None
    cc = None
    nu = None

    datafile = './fts_disk_center.idlsave'

    def __init__(self):
        
        # watt / (cm2 ster AA)  as emitted at solar surface
        t = readsav(self.datafile)
        
        # convert to J s-1 m-2 m-1 sr-1

        clight=2.99792458e8         #speed of light [m/s]                                  
        aa_to_m=1e-10                                                                        
        cm_to_m=1e-2                       
        
        self.ll = t['ftswav']* aa_to_m
        self.nu = clight / self.ll
        self.ii = t['ftsint'] * cm_to_m**(-2) * aa_to_m**(-1) # from from W /( cm2 ster AA) to  W / (m2 ster m)  
        self.cc = t['ftscnt'] * cm_to_m**(-2) * aa_to_m**(-1) 
       
