"""
Created on Wed Nov 27 12:42:17 2024

@author: esoria
@review: cruiz

This script is an example that reproduces the results for a zero order 6 layer quartz stack
"""
import numpy as np
from Retarder import Retarder
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

retarder = Retarder()

wl = np.arange(380, 1580, 2)  # wavelength range in nanometers

#Dictionary storing the Mueller Matrix elements
mueller_elements = {f"e{i}{j}": [] for i in range(1, 5) for j in range(1, 5)}
mueller_elements_DT = {f"e{i}{j}": [] for i in range(1, 5) for j in range(1, 5)}

for wave in wl:
    MMt = retarder.compute_mueller_stack(wave, bonding="contact")
    MMt_DT = retarder.compute_mueller_stack(wave, bonding="contact", DT=2.0)
    for i in range(4):
        for j in range(4):
            field_name = f"e{i+1}{j+1}"
            mueller_elements[field_name].append(MMt[i, j])
            mueller_elements_DT[field_name].append(MMt_DT[i, j])

fig, axs = plt.subplots(4, 4, figsize=(12, 10), constrained_layout=True)
fig.suptitle('Mueller matrix of a  zero-order P stack quartz retarder', fontsize=16)
for i in range(4):
    for j in range(4):
        ax = axs[i, j]
        field_name = f"e{i + 1}{j + 1}"
        ax.plot(wl, mueller_elements[field_name], linewidth=1.5)
        ax.plot(wl, mueller_elements_DT[field_name], linewidth=1.5)
        ax.set_title(f"M[{i + 1}, {j + 1}]")
        ax.set_xlabel('Wavelength (nm)')
        ax.set_ylabel('MM')

plt.show()