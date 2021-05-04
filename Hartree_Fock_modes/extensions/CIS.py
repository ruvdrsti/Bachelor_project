import psi4
import numpy as np
from scipy.linalg import eigh
numpy_memory = 4
psi4.set_memory(int(5e8))
class CISMolecule():
    def __init__(self, molecule):
        """
        Will set up some variables we will need from the Molecule object

        input
        molecule: a Molecule object from the compChem package
        constraint: True if you want a CUHF wavefunction
        """
        self.id = molecule
        self.occupied = self.id.alpha + self.id.beta
        self.available = self.id.integrals.nbf()*2
        self.virtual = self.available - self.occupied
        self.E_0 = molecule.E_0
    

    def getTwoElectronIntegrals(self, exchange=True):
        """
        returns two electron integrals in MO basis
        
        input:
        exchange: test parameter
        """
        # getting the two electron integrals in correct basis => we need it in MO basis
        tei = self.id.elrep # given in chemists notation

        #change the basis of the tei
        tei_int = np.kron(np.eye(2), tei)
        tei_big = np.kron(np.eye(2), tei_int.T)

        C = self.C
        tei_ao = tei_big.transpose(0, 2, 1, 3) - tei_big.transpose(0, 2, 3, 1) # accounts for both coulomb and exchange, switch to physisists notation 
        if exchange:
            tei_mo = np.einsum("pQRS,pP->PQRS", np.einsum("pqRS,qQ->pQRS", np.einsum("pqrS,rR->pqRS", np.einsum("pqrs,sS->pqrS", tei_ao, C, optimize=True), C, optimize=True), C, optimize=True), C, optimize=True)
        else:
            tei_mo = tei_ao
        return tei_mo
    

    def displayCISHamiltonian(self):
        """displays the CIS hamiltonian in MO basis"""
        # getting the orbital energies
        if self.id.mode == "rhf": 
            F_a = self.id.guessMatrix_a
            F_b = F_a
        else:
            F_a = self.id.guessMatrix_a
            F_b = self.id.guessMatrix_b
        F = np.block([[F_a, np.zeros(F_a.shape)], [np.zeros(F_b.shape), F_b]])
        overlap = np.kron(np.eye(2), self.id.overlap)
        epsilon, C = eigh(F, overlap)
        self.epsilon = epsilon
        self.C = C
        F = np.einsum("pq, qr, rs->ps", C.T, F, C, optimize=True)

        tei_mo = self.getTwoElectronIntegrals()
        
        #getting the excitations
        excitations = []
        for orbital in range(self.occupied): # for every occupied orbital
            for another_orbital in range(self.occupied, self.available): # we can make an excitation to every virtual orbital
                excitations.append((orbital, another_orbital))
        self.excitations = excitations
        # getting the hamiltonian
        dim = (self.occupied )*(self.virtual)
        H_cis = np.zeros((dim, dim))
        for row, excitation in enumerate(excitations):
            i, a = excitation
            for collumn, another_excitation in enumerate(excitations):
                j, b = another_excitation   
                H_cis[row, collumn] = self.id.getElectronicEnergy()*(i == j)*(a == b) + F[a, b]*(i==j) - F[i, j]*(a==b) + tei_mo[a, j, i, b] 
        
        H_cis = np.block([[self.id.getElectronicEnergy(), np.zeros((1,9))], [np.zeros((9,1)), H_cis]])
        return H_cis


    def CalculateExcitations(self):
        """setting up some properties needed for later"""
        ham = self.displayCISHamiltonian()
        
        self.ham = ham 
        excitation_energies, coefs = eigh(ham)
       
        self.coefs = coefs
        self.excitation_energies = excitation_energies


    def GetExitations(self, filepath="NoNameGiven", alternate=False):
        """Get the excitation energies and the contributions"""
        if filepath == "NoNameGiven":
            raise ValueError("no path specified")
        from pathlib import Path
        Path(f"{filepath}").touch()
        datafile = open(f"{filepath}", "w")
        self.CalculateExcitations()
        contrib = self.coefs**2
        energies = self.excitation_energies
        excitations = ["ground"] + self.excitations
        datafile.writelines(f"scf energy for {self.id.mode}: {self.id.getElectronicEnergy()}\n")
        for state, energy in enumerate(energies):
            datafile.writelines(f" {state + 1} : {energy}\n")
            for excitation, contribution in enumerate(contrib[:, state]):
                if contribution*100 > 1:
                    try:
                        datafile.writelines(f"\t{contribution:.3f}: {excitations[excitation]}\n")
                    except IndexError:
                        datafile.writelines(f"ground\n")
  
        datafile.close()
    
    
    def getCISEnergy(self, orbital=0):
        """
        will calculate the CIS energy
        
        input:
        orbital: you can enter the orbital from which you want the energy, default the lowest energy excitation
        """
        A = 0
        for number, excitation in enumerate(self.excitations):
            i, a = excitation
            A += (self.epsilon[a] - self.epsilon[i])*self.coefs[number, orbital]**2

        B = 0
        tei_mo = self.getTwoElectronIntegrals()
        for number, excitation in enumerate(self.excitations):
            i, a = excitation
            for another_number, another_excitation, in enumerate(self.excitations):
                j, b = another_excitation
                B += self.coefs[number, orbital]*self.coefs[another_number, orbital]*tei_mo[a, j, i, b]
        return A + B