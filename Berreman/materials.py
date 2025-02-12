# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 16:07:15 2024

@author: esoria
"""
import numpy as np

def fusedsilica(w):
    """
    Sellmeier dispersion equation for refractive index of fused silica.
    Parameters:
        w: Wavelength in nm (numpy array or scalar)
    Returns:
        Refractive index of fused silica (numpy array or scalar)
    """
    # Sellmeier coefficients
    B1 = 6.96166300E-01
    B2 = 4.07942600E-01
    B3 = 8.97479400E-01
    C1 = 4.67914826E-03
    C2 = 1.35120631E-02
    C3 = 9.79340025E+01
    # Sellmeier equation
    n = sellmeier(w, [B1, B2, B3, C1, C2, C3])
    return n

def quartz(w):
    """
    Laurent dispersion equation for crystalline quartz.

    Parameters:
        w: Wavelengths in nm (numpy array or scalar).

    Returns:
        A numpy array where each row corresponds to a wavelength and contains:
        [ne, no, d(ne)/dx, d(no)/dx].
    """
    w = w / 1000  # Convert wavelengths to micrometers
    ww = w ** 2

    # e-index coefficients
    B1_e = 2.38490E+00
    B2_e = -1.25900E-02
    B3_e = 1.07900E-02
    C1_e = 1.65180E-04
    C2_e = -1.94741E-06
    C3_e = 9.36476E-08
    ne = np.sqrt(
        B1_e + B2_e * ww + B3_e / ww + C1_e / ww ** 2 + C2_e / ww ** 3 + C3_e / ww ** 4
    )
    dnedx = (
            (B2_e * w - B3_e / w ** 3 - 2 * C1_e / w ** 5 - 3 * C2_e / w ** 7 - 4 * C3_e / w ** 9)
            / ne
            / 1000
    )

    # o-index coefficients
    B1_o = 2.35728E+00
    B2_o = -1.17000E-02
    B3_o = 1.05400E-02
    C1_o = 1.34143E-04
    C2_o = -4.45368E-07
    C3_o = 5.92362E-08
    no = np.sqrt(
        B1_o + B2_o * ww + B3_o / ww + C1_o / ww ** 2 + C2_o / ww ** 3 + C3_o / ww ** 4
    )
    dnodx = (
            (B2_o * w - B3_o / w ** 3 - 2 * C1_o / w ** 5 - 3 * C2_o / w ** 7 - 4 * C3_o / w ** 9)
            / no
            / 1000
    )

    # Combine results into a single numpy array
    return np.column_stack((ne, no, dnedx, dnodx))

def quartz_retarder(w):
    """
    index etc... for the quartz retarder
    :param w:
    :return:
    """
    #ordinary coefficients
    AO, BO, CO, DO, FO = 1.28604141, 1.07044083, 1.00585997e-2, 1.10202242, 100
    n_o = sellmeier(w, [AO, BO, CO, DO, FO])
    #extraordinary coefficients
    AE, BE, CE, DE, FE = 1.28851804, 1.09509924, 1.02101864e-2, 1.15662475, 100
    n_e = sellmeier(w, [AE, BE, CE, DE, FE])
    return n_o, n_e



def mgf2(w):
    """
    Sellmeier dispersion equation for crystalline magnesium fluoride (MgF2).

    Parameters:
        w: Wavelengths in nm (numpy array or scalar).

    Returns:
        A numpy array where each row corresponds to a wavelength and contains:
        [ne, no, d(ne)/dx, d(no)/dx].
    """
    w = w / 1000  # Convert wavelengths to micrometers
    ww = w ** 2

    # e-index coefficients
    B1_e = 4.13440230E-01
    B2_e = 5.04974990E-01
    B3_e = 2.49048620E+00
    C1_e = 1.35737865E-03
    C2_e = 8.23767167E-03
    C3_e = 5.65107755E+02
    ne = np.sqrt(
        1 + B1_e * ww / (ww - C1_e) + B2_e * ww / (ww - C2_e) + B3_e * ww / (ww - C3_e)
    )
    dnedx = -(
            B1_e * C1_e * w / (ww - C1_e) ** 2 +
            B2_e * C2_e * w / (ww - C2_e) ** 2 +
            B3_e * C3_e * w / (ww - C3_e) ** 2
    ) / ne / 1000

    # o-index coefficients
    B1_o = 4.87551080E-01
    B2_o = 3.98750310E-01
    B3_o = 2.31203530E+00
    C1_o = 1.88217800E-03
    C2_o = 8.95188847E-03
    C3_o = 5.66135591E+02
    no = np.sqrt(
        1 + B1_o * ww / (ww - C1_o) + B2_o * ww / (ww - C2_o) + B3_o * ww / (ww - C3_o)
    )
    dnodx = -(
            B1_o * C1_o * w / (ww - C1_o) ** 2 +
            B2_o * C2_o * w / (ww - C2_o) ** 2 +
            B3_o * C3_o * w / (ww - C3_o) ** 2
    ) / no / 1000

    # Combine results into a single numpy array
    return np.column_stack((ne, no, dnedx, dnodx))

def sellmeier(w, coeffs):
    """Computes the refraction index with the selmeier eq.s."""
    w /= 1000  # Converts to microns
    ww = w**2
    n = coeffs[0]+(coeffs[1] * ww)/(ww - coeffs[2])+(coeffs[3] * ww)/(ww-coeffs[4])
    return np.sqrt(n)


def birrefringence(w, material):
    An = material["H"] + (material["I"] * w**2)/(w**2 - material["G"]) + (material["J"] * w**2)/(w**2 - material["L"])
    return An