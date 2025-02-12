"""
Created on Wed Nov 27 12:42:17 2024

@author: esoria
@review: cruiz

This script is an example that reproduces the results for a zero order 6 layer quartz stack
"""
import numpy as np
from Retarder import Retarder
import os
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import pandas as pd

# Create directories for saving results
os.makedirs("plots", exist_ok=True)
os.makedirs("data", exist_ok=True)

retarder = Retarder()
R=150000 #max resolving power

pol_waves = [393.3, 396.8, 422.7, 453, 455.4, 460.7, 492.3, 514, 517.3,
                     525, 587.6, 589, 617.3, 630, 656.3, 709, 815.1, 849.8, 854,
                     1083.0, 1565.0]

for pol in pol_waves:
    dl = pol/(R*2)
    wl = np.arange(pol-1, pol+1, dl)  # +/- 1 nm interval around the target line

    # Dictionary storing the Mueller Matrix elements
    mueller_elements = {f"e{i}{j}": [] for i in range(1, 5) for j in range(1, 5)}
    # Dictionary storing the Mueller Matrix elements
    mueller_elements_DT = {f"e{i}{j}": [] for i in range(1, 5) for j in range(1, 5)}

    for wave in wl:
        MMt = retarder.compute_mueller_stack(wave, bonding="contact")
        MMt_DT = retarder.compute_mueller_stack(wave, bonding="contact", DT=2.0)
        for i in range(4):
            for j in range(4):
                field_name = f"e{i+1}{j+1}"
                mueller_elements[field_name].append(MMt[i, j])
                mueller_elements_DT[field_name].append(MMt_DT[i, j])

        # Save data to CSV
    df = pd.DataFrame({"Wavelength": wl})
    for key in mueller_elements:
        df[key] = mueller_elements[key]
        df[f"{key}_DT"] = mueller_elements_DT[key]

    csv_filename = f"data/mueller_matrix_{pol}.csv"
    df.to_csv(csv_filename, index=False)
    print(f"Data saved: {csv_filename}")

    fig, axs = plt.subplots(4, 4, figsize=(12, 10), constrained_layout=True)
    fig.suptitle('Mueller matrix of a  zero-order P stack quartz retarder ' +str(pol), fontsize=16)
    for i in range(4):
        for j in range(4):
            ax = axs[i, j]
            field_name = f"e{i + 1}{j + 1}"
            ax.plot(wl, mueller_elements[field_name], linewidth=1.5, label="No DT")
            ax.plot(wl, mueller_elements_DT[field_name], linewidth=1.5, label="DT=2.0")
            ax.set_title(f"M[{i + 1}, {j + 1}]")
            ax.set_xlabel('Wavelength (nm)')
            ax.set_ylabel('MM')
            ax.legend()


    # Save the figure
    plot_filename = f"plots/mueller_matrix_{pol}.png"
    plt.savefig(plot_filename, dpi=300)
    plt.close()
    print(f"Plot saved: {plot_filename}")
    del(mueller_elements)
    del(mueller_elements_DT)