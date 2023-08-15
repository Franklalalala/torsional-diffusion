import os

import ase
import numpy as np
import openpyxl
import rdkit
from ase.atom import Atom
from ase.atoms import Atoms
from ase.db import connect
from ase.io import read, write
from ase.neighborlist import build_neighbor_list
from ase.visualize import view
from openpyxl import load_workbook
from rdkit import Chem
from rdkit import Geometry
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdmolfiles
from rdkit.Chem import rdmolops


def atom_2_mol(an_atoms: ase.atoms.Atoms):
    write(filename='temp.xyz', images=an_atoms)
    a_mol = rdmolfiles.MolFromXYZFile('temp.xyz')
    os.remove('temp.xyz')
    return a_mol


def set_mol_prop_from_row(row, mol):
    mol.SetProp('ep_id', f'{str(row.ep_id)}')
    mol.SetProp('smile', f'{row.smile}')
    mol.SetProp('binding_e', f'{str(row.binding_e)}')
    mol.SetProp('viscosity', f'{row.viscosity}')
    mol.SetProp('dielectric_constant', f'{row.dielectric_constant}')
    mol.SetProp('lumo', f'{row.lumo}')
    mol.SetProp('homo', f'{row.homo}')
    return mol


def get_atoms_from_mol(mol: rdkit.Chem.rdchem.Mol):
    conf = mol.GetConformer()
    an_atoms = Atoms()
    for i in range(conf.GetNumAtoms()):
        position = conf.GetAtomPosition(i)
        atom = mol.GetAtoms()[i]
        a_symbol = atom.GetSymbol()
        an_new_atom = Atom(symbol=a_symbol, position=(position.x, position.y, position.z))
        an_atoms.append(an_new_atom)
    return an_atoms


def edit_ase_db(ase_db_path):
    with connect(ase_db_path) as db:
        meta = db.metadata
        if not "_distance_unit" in meta.keys():
            meta["_distance_unit"] = 'Ang'
        if not "_property_unit_dict" in meta.keys():
            meta["_property_unit_dict"] = {
                'viscosity': "Hartree",
                'dielectric_constant': "Hartree",  # relative
                'binding_e': "eV",
                'homo': "eV",
                'lumo': "eV",
            }
        if "atomrefs" not in meta.keys():
            meta["atomrefs"] = {}
        db.metadata = meta


def xlsx2db(db_path: str, xlsx_path: str, sheet_name: str):
    wb = load_workbook(filename=xlsx_path)
    ws = wb[sheet_name]

    real_count = 0
    fail_count = 0
    with connect(db_path) as db:
        for i, row in enumerate(ws.values):
            if i == 0:
                continue
            # if i>10:
            #     continue
            a_smile = row[2]
            try:
                a_mol = Chem.MolFromSmiles(a_smile)
                a_mol_with_H = Chem.AddHs(a_mol)
                AllChem.EmbedMolecule(a_mol_with_H, useRandomCoords=True, maxAttempts=1000000)
                AllChem.MMFFOptimizeMolecule(a_mol_with_H)

                an_atoms = get_atoms_from_mol(mol=a_mol_with_H)

                an_id = row[1]
                an_id = int(an_id.split('-')[1])
                a_binding_e = float(row[3])
                a_dielectric_c = np.log10(float(row[4]))
                a_viscosity = np.log10(float(row[5]))
                a_lumo = float(row[6])
                a_homo = float(row[7])

                db.write(atoms=an_atoms, binding_e=a_binding_e, viscosity=a_viscosity,
                         dielectric_constant=a_dielectric_c, lumo=a_lumo, homo=a_homo,
                         smile=a_smile, ep_id=an_id)
                real_count = real_count + 1
                # print(f'real: {real_count}')
            except Exception as e:
                print(e)
                fail_count = fail_count + 1
                print(f'real: {real_count}')
                print(f'fail: {fail_count}')
                with open('fail_smile', 'a') as f_f:
                    f_f.write(a_smile)
                    f_f.write('\n')
    edit_ase_db(db_path)

# xlsx2db(db_path=r'iid_test.db', xlsx_path=r'1_CHO_47371_uninf_20230706_iid_test.xlsx', sheet_name='Sheet1')




