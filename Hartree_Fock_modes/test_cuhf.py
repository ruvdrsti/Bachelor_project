from cuhf import CUHFMolecule
from uhf import UHFMolecule
import psi4
import numpy as np
psi4.set_options({"basis":"sto-3g"})
distances = np.arange(0.2, 5, 0.1)
energies_cuhf = []
spincont_cuhf = []
for distance in distances:
    h2 = CUHFMolecule(f"""
    H 0 0 0 
    H 0 0 {distance}
    """)
    h2.setConvergence(1e-12)
    E = h2.iterator(mute=True, criterion="energy", mixedGuess=True)
    energies_cuhf.append(E[0])
    spincont_cuhf.append(h2.getSpinContamination())

print(h2.displayFockMatrix("alpha"))
import matplotlib.pyplot as plt
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(distances, energies_cuhf, label="cuhf energy")
ax2.scatter(distances, spincont_cuhf, label="spin contamination", color="black", marker=".")
ax1.set_xlabel("distance in Angstrom")
ax1.set_ylabel("energy in Hartree")
ax2.set_ylabel("spin contamination")
ax1.hlines(0, 0, 5, color="red", linestyles="--")
fig.legend( loc="upper right")
ax1.axis([0, 5, -1.3, 1.1])
ax2.axis([0, 5, -1.3, 1.1])
plt.savefig("test")
