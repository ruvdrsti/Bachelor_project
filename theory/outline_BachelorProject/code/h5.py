import psi4
import numpy as np
from compChem import Hartree_fock
import sys
sys.path.insert(0, "./../../../Hartree_Fock_modes/extensions/")
sys.path.insert(0, "./../../../Hartree_Fock_modes/")
from CIS import CISMolecule
from cuhf import CUHFMolecule
from uhf import UHFMolecule
psi4.set_options({"basis":"sto-3g", "reference":"cuhf", "scf_type":"pk"})
e_list = []
distances = np.arange(0.2, 3, 0.1)
distance = 1
geometry = f"""
H 0 0 0 
H 1 0 0
H {1+np.cos(1.2566370614)} {np.sin(1.2566370614)} 0
H {2*np.cos(0.62831853072)*np.cos(1.2566370614)} {2*np.cos(0.62831853072)*np.sin(1.2566370614)} 0
H {-np.sin(0.31415926536)} {np.cos(0.31415926536)} 0
"""
h5 = UHFMolecule(geometry)
h5.setConvergence(1e-12)
E = h5.iterator(mute=True, criterion="energy", mixedGuess=True)
e_list.append(E[0])
cis = CISMolecule(h5)
cis.GetExitations("h5_4.txt")

