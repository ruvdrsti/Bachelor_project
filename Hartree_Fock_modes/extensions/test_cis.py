from compChem.Hartree_Fock_modes.cuhf import CUHFMolecule
from CIS import CISMolecule
import psi4
# molecule n°2: allyl radical
psi4.set_options({"basis":"sto-3g", "scf_type":"pk", "reference":"cuhf", "e_convergence":"1e-12"})
allyl = CUHFMolecule("""
0 1
O
H 1 1.1
H 1 1.1 2 104
symmetry c1
""")
allyl.setConvergence(1e-12)
end_data = allyl.iterator(mute=True, criterion="energy", mixedGuess=True)
print(f"scf energy {end_data[0]: .14f}, {end_data[1]} iterations")
cis = CISMolecule(allyl)
cis.GetExitations("extensions/allyl_test.txt")
cis.plotExcitations("extensions/test_cis.png")