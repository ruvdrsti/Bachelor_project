import psi4

import numpy as np
from cuhf import CUHFMolecule
from uhf import UHFMolecule
psi4.core.set_output_file('CUHF.dat', False)
psi4.set_memory("3 GB")
psi4.set_options({"basis":"sto-3g", "scf_type":"pk", "reference":"cuhf", "cubeprop_tasks":"orbitals", "e_convergence":"1e-12"})
h2o = CUHFMolecule("""
H 0 0 0
H 0 0.86602540378 0.5
H 0 0 1
units angstrom""")
h2o.setConvergence(1e-12)
E, wfn_RHF = psi4.energy("scf", mol=h2o.id, return_wfn=True)
F = h2o.iterator(mute=True, criterion="energy")
Ca = h2o.getEigenStuff("alpha")[1]
Cb = h2o.getEigenStuff("beta")[1]
print(Ca)
print(Cb)
#Ca[:,0] = 0.1
#Cb[:,0] = 0.1
wfn_RHF.Ca().copy(psi4.core.Matrix.from_array(Ca))
wfn_RHF.Cb().copy(psi4.core.Matrix.from_array(Cb))
psi4.cubeprop(wfn_RHF)





