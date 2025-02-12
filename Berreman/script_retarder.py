# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:42:17 2024

@author: esoria

This script is an example that reproduces the results of the figure Fig. 5 Harrington et al.: Polarization modeling and
predictions for DKIST part 2: application.
"""
import numpy as np
from Berreman import Berreman
import matplotlib.pyplot as plt
import materials as mat
#import Berreman.materials as mat

pi = np.pi

B = Berreman()

theta_i = 0  # Incidence angle

wl = np.arange(380, 900, 2)  # wavelength range in nanometers

#Dictionary storing the Mueller Matrix elements
mueller_elements = {f"e{i}{j}": [] for i in range(1, 5) for j in range(1, 5)}

materials = {
    "Quartz": {"H": 0.61339 * 10 ** -3, "I": 8.2187 * 10 ** -3, "G": 13.476 * 10 ** -3, "J": 13.644 * 10 ** -3,
               "L": 82.208},
    "Sapphire": {"H": 1.8299 * 10 ** -3, "I": -9.6436 * 10 ** -3, "G": 10.432 * 10 ** -3, "J": -0.38758 * 10 ** -3,
                 "L": 21.660},
    "MgF2": {"H": -19.364 * 10 ** -3, "I": 30.992 * 10 ** -3, "G": 2.3253 * 10 ** -3, "J": 40.060 * 10 ** -3,
             "L": 388.37}
    }

# air + MgF2 + Quartz + air
for wave in wl:
    wave /= 1000 #converts to microns
    no = -0.4041 * wave ** 5 + 1.5087 * wave ** 4 - 2.2634 * wave ** 3 + 1.7245 * wave ** 2 - 0.6841 * wave + 1.492
    b = mat.mgf2(wave*1000)
    beta = 1 * np.sin(np.radians(theta_i))
    An = mat.birrefringence(wave, materials["MgF2"])
    thickC = 13.587 * 1000 / (wave * 1000)  # Thickness (109nm)

    Fc = B.fmat([1, 1, 1, 0, 0, np.inf], beta)
    Fs = B.fmat([1, 1, 1, 0, 0, np.inf], beta)
    M = B.cmat([no, no + An, no, 0, 0, 0, thickC], beta)

    MMr, MMt, Sr, St = B.cmmat(Fc, M, Fs)
    norm = MMt[0, 0]
    MMt = MMt / norm
    MMt[0, 0] = norm

    # Stores the elements of the Mueller Matrix
    for i in range(4):
        for j in range(4):
            field_name = f"e{i+1}{j+1}"
            mueller_elements[field_name].append(MMt[i, j])

# Crates 4x4 figure template
fig, axs = plt.subplots(4, 4, figsize=(12, 10), constrained_layout=True)
fig.suptitle('Mueller matrix of a  zero-order quartz retarder', fontsize=16)
# Plots the Matrix elements
for i in range(4):
    for j in range(4):
        ax = axs[i, j]
        field_name = f"e{i + 1}{j + 1}"
        ax.plot(wl, mueller_elements[field_name], linewidth=1.5)
        ax.set_title(f"M[{i + 1}, {j + 1}]")
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('MM')

            # Configures the specific limits for the graphs

        if j == 0 and i == 0:
            ax.set_ylim([0.86, 1])
        elif j == 0 and i == 1:
            ax.set_ylim([-0.1, 0.1])
        elif j == 1 and i == 0:
            ax.set_ylim([-0.1, 0.1])
        elif j == 2 and i == 3:
            ax.set_ylim([0.7, 1])

        elif j == 3 and i == 2:
            ax.set_ylim([-1, -0.7])
        elif j == 2 and i == 2:
            ax.set_ylim([-0.6, 0.6])
        elif j == 3 and i == 3:
            ax.set_ylim([-0.6, 0.6])
        else:
            ax.set_ylim([-1, 1])

        ax.set_xlim([380, 900])
        ax.grid(True)

plt.show()