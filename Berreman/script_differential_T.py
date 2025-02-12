"""


This script is an example that reproduces the results for a zero order 6 layer quartz stack and computes the effects
of differential temperature variations
"""
import numpy as np
from Retarder import Retarder
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

R=150000 #max resolving power

pol_waves = [393.3, 396.8, 422.7, 453, 455.4, 460.7, 492.3, 514, 517.3,
                     525, 587.6, 589, 617.3, 630, 656.3, 709, 815.1, 849.8, 854,
                     1083.0, 1565.0]

retarder = Retarder()



for pol in pol_waves:
    dl = pol/(R*2)
    wl = np.arange(pol-1,pol+1, dl)  # +/- 1 nm interval around the target line

    # Dictionary storing the Mueller Matrix elements
    mueller_elements = {f"e{i}{j}": [] for i in range(1, 5) for j in range(1, 5)}
    # Dictionary storing the Mueller Matrix elements
    mueller_elements_DT = {f"e{i}{j}": [] for i in range(1, 5) for j in range(1, 5)}

    for wave in wl:
        MMt = retarder.compute_mueller_stack(wave, bonding="contact")
        MMt_DT = retarder.compute_mueller_stack(wave, bonding="contact", DT=2.0)
        for i in range(4):
            for j in range(4):
                field_name = f"e{i + 1}{j + 1}"
                mueller_elements[field_name].append(MMt[i, j])
                mueller_elements_DT[field_name].append(MMt_DT[i, j])

    fig, axs = plt.subplots(4, 4, figsize=(12, 10), constrained_layout=True)
    fig.suptitle('Mueller matrix of a  zero-order P stack quartz retarder'+str(pol), fontsize=16)
    for i in range(4):
        for j in range(4):
            ax = axs[i, j]
            field_name = f"e{i + 1}{j + 1}"
            ax.plot(wl, np.asarray(mueller_elements[field_name])-np.asarray(mueller_elements_DT[field_name]), linewidth=1.5)
            ax.set_title(f"M[{i + 1}, {j + 1}]")
            ax.set_xlabel('Wavelength (nm)')
            ax.set_ylabel('MM')
    plt.show()
    del(mueller_elements)
    del(mueller_elements_DT)