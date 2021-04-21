from hf_backbone import Molecule
from cuhf import CUHFMolecule
import psi4

# molecule n°2: allyl radical
psi4.set_options({"basis":"sto-3g", "scf_type":"pk", "reference":"cuhf", "e_convergence":"1e-12"})
allyl = CUHFMolecule("""
0 2
H
C 1 r2
C 2 r3 1 a3
C 2 r3 1 a3 3 180.
H 3 r5 2 a5 1 0.
H 4 r5 2 a5 1 0.
H 3 r7 2 a7 1 180.
H 4 r7 2 a7 1 180.

r2=1.08424658
r3=1.40526604
r5=1.08095381
r7=1.08131649
a3=117.99450641
a5=121.41544408
a7=121.21891262
symmetry c1

""")
allyl.setConvergence(1e-12)
end_data = allyl.iterator(mute=True, criterion="energy", mixedGuess=True)
print(f"scf energy {end_data[0]: .14f}, {end_data[1]} iterations")