class CISMolecule(CISMolecule):
    def alternativeHamiltonian(self):
        """alternative approach to the hamiltonian"""
        # getting the orbital energies
        epsilon_a, C_a = eigh(self.id.displayFockMatrix("alpha"), self.id.overlap)
        epsilon_b, C_b = eigh(self.id.displayFockMatrix("beta"), self.id.overlap)
        epsilon = np.append(epsilon_a, epsilon_b)
        sortedorder = np.argsort(epsilon)
        epsilon.sort()
        self.epsilon = epsilon
        
        # make the C matrix => it contains all the orbitals, both alpha and beta
        C = np.block([[C_a, np.zeros(C_a.shape)], [np.zeros(C_b.shape), C_b]])
        C = C[:, sortedorder] # sort the eigenfunctions
        self.C = C

        tei_mo = self.getTwoElectronIntegrals()
        tei_mo_no_exch = self.getTwoElectronIntegrals(exchange=False)

        # get excitations alpha -> alpha
        excitations_aa = []
        for orbital in range(self.id.alpha):
            for another_orbtial in range(self.id.alpha, self.id.integrals.nbf()):
                 excitations_aa.append((orbital, another_orbtial))

        # get excitations alpha -> beta 
        excitations_ab = []
        for orbital in range(self.id.alpha):
            for another_orbtial in range(self.id.beta, self.id.integrals.nbf()):
                 excitations_ab.append((orbital, another_orbtial))

        # get excitations beta -> alpha
        excitations_ba = []
        for orbital in range(self.id.beta):
            for another_orbtial in range(self.id.alpha, self.id.integrals.nbf()):
                 excitations_aa.append((orbital, another_orbtial))    

        # get excitations beta -> beta
        excitations_bb = []
        for orbital in range(self.id.beta):
            for another_orbtial in range(self.id.beta, self.id.integrals.nbf()):
                 excitations_aa.append((orbital, another_orbtial))
        
        HCis = np.zeros((self.occupied*self.virtual, self.occupied*self.virtual))
        add = self.id.integrals.nbf()
        # generate aa block
        for row, excitation in enumerate(excitations_aa):
            i, a = excitation
            for collumn, another_excitation in enumerate(excitations_aa):
                j, b = another_excitation
                HCis[row, collumn] = (epsilon_a[a] - epsilon_a[i])*(i == j)*(a == b)+ tei_mo[a, j, i, b]
        
        # generate bb block
        for row, excitation in enumerate(excitations_bb):
            i, a = excitation
            for collumn, another_excitation in enumerate(excitations_bb):
                j, b = another_excitation
                HCis[row + add, collumn + add] = (epsilon_b[a] - epsilon_b[i])*(i == j)*(a == b) + tei_mo[a + add, j + add, i + add, b + add]
        
        # generate ab block
        for row, excitation in enumerate(excitations_ab):
            i, a = excitation
            for collumn, another_excitation in enumerate(excitations_ab):
                j, b = another_excitation
                HCis[row, collumn + add] = (epsilon_b[a] - epsilon_a[i])*(i == j)*(a == b) + tei_mo[a + add, j, i, b + add]
        
        # generate ba block
        for row, excitation in enumerate(excitations_ba):
            i, a = excitation
            for collumn, another_excitation in enumerate(excitations_ba):
                j, b = another_excitation
                HCis[row + add, collumn] = (epsilon_a[a] - epsilon_b[i])*(i == j)*(a == b) + tei_mo[a, j + add, i + add, b]
        
        return HCis