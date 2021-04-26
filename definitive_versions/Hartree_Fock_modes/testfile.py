from hf_backbone import Molecule
from cuhf import CUHFMolecule
import psi4

psi4.set_options({"basis":"sto-3g", "scf_type":"pk", "reference":"cuhf", "e_convergence":"1e-12"})
hydrogen = CUHFMolecule("""H""")
hydrogen.setConvergence(1e-12)
h_energy = hydrogen.iterator(mute=True)

import numpy as np
distances = np.arange(0.2, 5, 0.1)
energies = []
for distance in distances:
    h2 = CUHFMolecule(f"""
    H 0 0 0 
    H 0 0 {distance}
    """)
    h2.setConvergence(1e-12)
    E = h2.iterator(mute=True, criterion="energy", mixedGuess=True)
    energies.append(E[0] - h_energy[0]*2)

import matplotlib.pyplot as plt
plt.plot(distances, energies, label="cuhf energy")
plt.hlines(0, 0, 5, color="red", label="zero")
plt.xlabel("distance in Angstrom")
plt.ylabel("energy in Hartree")
plt.legend(loc="upper right")
plt.savefig("temp.png")