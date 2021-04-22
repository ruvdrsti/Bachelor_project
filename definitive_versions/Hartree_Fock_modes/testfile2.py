from rhf import RHFMolecule
from extensions.CIS import CISMolecule
import psi4

psi4.set_options({"basis":"sto-3g", "scf_type":"pk", "reference":"rhf", "e_convergence":"1e-12"})
allyl = RHFMolecule("""
0 1
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")
allyl.setConvergence(1e-12)
end_data = allyl.iterator(mute=True, criterion="energy")
cis = CISMolecule(allyl)
cis.GetExitations("testdata.txt")