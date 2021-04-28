# this file is meant to demonstrate the different data between CIS measurements.
import psi4
from compChem import Hartree_fock
from compChem.CIS import CISMolecule
test = open(r"/user/gent/440/vsc44013/Bachelor_project/Bachelor_project/theory/outline_BachelorProject/data/test_cis_data/guessmatrices_force.txt", "w")
psi4.set_options({"basis":"sto-3g", "reference":"rhf"})
water = Hartree_fock.RHFMolecule("""
0 1
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")
water.setConvergence(1e-12)
end_data = water.iterator(criterion="energy", mute=True)
cis = CISMolecule(water)
print(cis.id.guessMatrix_a)
end_data = water.iterator(criterion="energy", mute=True)
cis = CISMolecule(water)
print(cis.id.guessMatrix_a)