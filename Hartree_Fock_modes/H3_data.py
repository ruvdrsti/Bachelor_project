import psi4

import numpy as np

psi4.core.set_output_file('output.dat', False)

psi4.set_options({"basis":"sto-3g", "scf_type":"pk", "reference":"rohf", "ex_level":1, "num_roots":2})
h2o = psi4.geometry("""
H 0 0 0
H 0 0.86602540378 0.5
H 0 0 1
units angstrom""")
psi4.energy("DETCI")



