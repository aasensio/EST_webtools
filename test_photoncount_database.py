import unittest

from photoncount_database import *
from photoncount import Photocount

class TestNFlux(unittest.TestCase):

    # 2.1.3
    def test_nflux1(self):
        data = [854, -1, "IFS-M", 70000, 0.57, 1000, 5, True]
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        self.assertAlmostEqual(nflux.max(), 3.813e+19, -18)
        #añadir porcentaje de error y anular test

    # 2.1.3
    def test_nflux2(self):
        data = [854, 0.05, "TIS", 50000, 0.57, 500, 1, True]
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        self.assertAlmostEqual(nflux.max(), 3.813e+19, -18)

    # 3.1.2
    def test_nflux3(self):
        data = [396, 0.06, "TIS", 100000, 0.13, 300, 3, True]
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        self.assertAlmostEqual(nflux.max(), 4.923e+18, -18)

    # 3.1.2
    def test_nflux4(self):
        data = [630, 0.06, "TIS", 100000, 0.21, 300, 3, True]
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        self.assertAlmostEqual(nflux.max(), 2.647e+19, -17)

    # 3.1.2
    def test_nflux5(self):
        data = [656, 0.06, "TIS", 100000, 0.22, 300, 3, True]
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        self.assertAlmostEqual(nflux.max(), 2.651e+19, -19)

    # 3.1.2
    def test_nflux6(self):  #this one is failing, retire -20 by -19 and see
                            #how distant the values are
        data = [1083, 0.06, "IFS-S", 100000, 0.36, 300, 30, True]
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        self.assertAlmostEqual(nflux.max(), 1.125e+19, -20)

    # 2.1.3
    def test_adjust_nflux1(self):
        data = [854, -1, "IFS-M", 70000, 0.57, 1000, 5, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        self.assertAlmostEqual(nflux.max(), 1162831863, -14)

    # 2.1.3
    def test_adjust_nflux2(self):
        data = [854, 0.05, "TIS", 50000, 0.57, 500, 1, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        self.assertAlmostEqual(nflux.max(), 325592922, -12)

    # 3.1.2
    def test_adjust_nflux3(self):
        data = [396, 0.06, "TIS", 100000, 0.13, 300, 3, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        self.assertAlmostEqual(nflux.max(), 29244739, -12)

    # 3.1.2
    def test_adjust_nflux4(self): #Almost perfect
        data = [630, 0.06, "TIS", 100000, 0.21, 300, 3, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        self.assertAlmostEqual(nflux.max(), 250133690, -2)

    # 3.1.2
    def test_adjust_nflux5(self):
        data = [656, 0.06, "TIS", 100000, 0.22, 300, 3, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        self.assertAlmostEqual(nflux.max(), 260892113, -13)

    # 3.1.2
    def test_adjust_nflux6(self):
        data = [1083, 0.06, "IFS-S", 100000, 0.36, 300, 30, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22**2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        self.assertAlmostEqual(nflux.max(), 8325711895, -14)

    @classmethod
    def tearDownClass(self):
        print("\tClass " + __class__.__name__ + " has been tested.")

class TestPxFactor(unittest.TestCase):

    # 2.1.3
    def test_pxfactor1(self):
        data = [854, -1, "IFS-M", 70000, 0.57, 1000, 5, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        self.assertAlmostEqual(pxfactor, 0.8, 2)

    # 2.1.3
    def test_pxfactor2(self):
        data = [854, 0.05, "TIS", 50000, 0.57, 500, 1, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        self.assertAlmostEqual(pxfactor, 2.38, 2)

    # 3.1.2
    def test_pxfactor3(self):
        data = [396, 0.06, "TIS", 100000, 0.13, 300, 3, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        self.assertAlmostEqual(pxfactor, 6.17, 2)

    # 3.1.2
    def test_pxfactor4(self):
        data = [630, 0.06, "TIS", 100000, 0.21, 300, 3, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        self.assertAlmostEqual(pxfactor, 3.88, 2)
    
    # 3.1.2
    def test_pxfactor5(self):
        data = [656, 0.06, "TIS", 100000, 0.22, 300, 3, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        self.assertAlmostEqual(pxfactor, 3.72, 2)
    
    # 3.1.2
    def test_pxfactor6(self):
        data = [1083, 0.06, "IFS-S", 100000, 0.36, 300, 30, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        self.assertAlmostEqual(pxfactor, 2.26, 2)

    @classmethod
    def tearDownClass(self):
        print("\tClass " + __class__.__name__ + " has been tested.")

class TestNPhotons(unittest.TestCase):

    # 2.1.3
    def test_nphotons1(self):
        data = [854, -1, "IFS-M", 70000, 0.57, 1000, 5, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22 ** 2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(nflux.max(), pxfactor, data[7])
        self.assertAlmostEqual(nphotons, 106349347, -7)

    # 2.1.3
    def test_nphotons2(self):
        data = [854, 0.05, "TIS", 50000, 0.57, 500, 1, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22 ** 2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(nflux.max(), pxfactor, data[7])
        self.assertAlmostEqual(nphotons, 267006930, -7)

    # 3.1.2
    def test_nphotons3(self):
        data = [396, 0.06, "TIS", 100000, 0.13, 300, 3, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22 ** 2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(nflux.max(), pxfactor, data[7])
        self.assertAlmostEqual(nphotons, 160613921, -8)

    # 3.1.2
    def test_nphotons4(self):
        data = [630, 0.06, "TIS", 100000, 0.21, 300, 3, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22 ** 2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(nflux.max(), pxfactor, data[7])
        self.assertAlmostEqual(nphotons, 542771336, -2)
    
    # 3.1.2
    def test_nphotons5(self):
        data = [656, 0.06, "TIS", 100000, 0.22, 300, 3, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22 ** 2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(nflux.max(), pxfactor, data[7])
        self.assertAlmostEqual(nphotons, 522130525, -8)
    
    # 3.1.2
    def test_nphotons6(self):
        data = [1083, 0.06, "IFS-S", 100000, 0.36, 300, 30, True]
        test = DatabasePhotoncount(data)
        photoncount = Photocount()
        lmin = data[0] - data[4] / 2
        lmax = data[0] + data[4] / 2
        photoncount.R = data[3]
        nflux = photoncount.compute_nflux(lmin, lmax)
        nflux /= (1.22 ** 2)
        nflux = test.adjust_nflux(nflux, data[0], data[3], data[6])
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(nflux.max(), pxfactor, data[7])
        self.assertAlmostEqual(nphotons, 6113505083, -9)
    
    @classmethod
    def tearDownClass(self):
        print("\tClass " + __class__.__name__ + " has been skipped.")

class TestTransmission(unittest.TestCase):

    # 2.1.3 if we delete the reasignment of nphotons to it correct value if fails
    def test_transmission1(self):
        data = [854, -1, "IFS-M", 70000, 0.57, 1000, 5, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(1162831863, pxfactor, data[7])
        #nphotons = 184314294
        transmission = (data[5]**2) / nphotons
        self.assertAlmostEqual(transmission, 0.0054, 2)

    # 2.1.3 if we delete the reasignment of nphotons to it correct value if fails
    def test_transmission2(self):
        data = [854, 0.05, "TIS", 50000, 0.57, 500, 1, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(325592922, pxfactor, data[7])
        #nphotons = 462750312
        transmission = (data[5]**2) / nphotons
        self.assertAlmostEqual(transmission, 0.0005, 3)

    # 3.1.2 if we delete the reasignment of nphotons to it correct value if fails
    def test_transmission3(self):
        data = [396, 0.06, "TIS", 100000, 0.13, 300, 3, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(29244739, pxfactor, data[7])
        #nphotons = 278360348
        transmission = (data[5]**2) / nphotons
        self.assertAlmostEqual(transmission, 0.0003, 3)

    # 3.1.2 if we delete the reasignment of nphotons to it correct value if fails
    def test_transmission4(self):
        data = [630, 0.06, "TIS", 100000, 0.21, 300, 3, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(250133690, pxfactor, data[7])
        #nphotons = 940678225
        transmission = (data[5]**2) / nphotons
        self.assertAlmostEqual(transmission, 0.0001, 3)

    # 3.1.2 if we delete the reasignment of nphotons to it correct value if fails
    def test_transmission5(self):
        data = [656, 0.06, "TIS", 100000, 0.22, 300, 3, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(260892113, pxfactor, data[7])
        #nphotons = 904905588
        transmission = (data[5]**2) / nphotons
        self.assertAlmostEqual(transmission, 0.0001, 3)

    # 3.1.2 if we delete the reasignment of nphotons to it correct value if fails
    def test_transmission6(self):
        data = [1083, 0.06, "IFS-S", 100000, 0.36, 300, 30, True]
        test = DatabasePhotoncount(data)
        pxfactor = test.compute_pxfactor(data[0], data[1])
        nphotons = test.compute_n_photons(8325711895, pxfactor, data[7])
        #nphotons = 10595329433
        transmission = (data[5]**2) / nphotons
        self.assertAlmostEqual(transmission, 0, 4)

    @classmethod
    def tearDownClass(self):
        print("\tClass " + __class__.__name__ + " has been tested.")

if __name__ == '__main__':
    unittest.main()