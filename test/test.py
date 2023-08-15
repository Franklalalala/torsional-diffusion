import sys

sys.path.append(r'F:\test_torsional_diffusion\torsional-diffusion\utils')

from ase_tools import *
from ase.db import connect

with connect('iid_test.db') as db:
    for a_row in db.select():
        an_atom = db.get_atoms(ep_id=a_row.ep_id)
        a_mol = atom_2_mol(an_atoms=an_atom)
        set_mol_prop_from_row(row=a_row, mol=a_mol)
        a=1