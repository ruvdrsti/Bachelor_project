from rhf import RHFMolecule
from uhf import UHFMolecule
from cuhf import CUHFMolecule
from extensions.CIS import CISMolecule
import psi4
import numpy as np

psi4.set_options({"basis":"sto-3g", "scf_type":"pk", "reference":"cuhf", "e_convergence":"1e-12"})
allyl = CUHFMolecule("""
H 0 0 0
H 0 0.86602540378 0.5
H 0 0 1
units angstrom""")
allyl.setConvergence(1e-12)
end_data = allyl.iterator(mute=True, criterion="energy")
cis = CISMolecule(allyl)
cis.GetExitations("testdata.txt")