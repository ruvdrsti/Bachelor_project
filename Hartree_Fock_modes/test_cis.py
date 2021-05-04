import psi4
from compChem import Hartree_fock
from cuhf import CUHFMolecule
from uhf import UHFMolecule
from extensions.CIS import CISMolecule
psi4.set_options({"basis":"sto-3g", "reference":"rohf", "scf_type":"pk", "num_roots":2, "ex_level":1, "diag_method":"davidson"})
psi4.core.set_output_file('ROHF.dat', False)
h3 = CUHFMolecule("""
H 0 0 0
H 0 0.86602540378 0.5
H 0 0 1
units angstrom""")
h3.setConvergence(1e-12)

E = h3.iterator(mute=True, criterion="energy")
cis = CISMolecule(h3)
cis.GetExitations("test_cis.txt")

