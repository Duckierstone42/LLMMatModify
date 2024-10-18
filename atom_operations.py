from ase import Atoms
import os
from ase.data import chemical_symbols
from ase.build import molecule
from pymatgen.io.ase import AseAtomsAdaptor
import numpy as np
from dotenv import load_dotenv
from mp_api.client import MPRester
import random

load_dotenv()

def convert_to_atom(chemical_formula):
    with MPRester(os.getenv("MP_API_KEY")) as mpr:
            docs = mpr.materials.summary.search(
                formula=chemical_formula,
                fields=["structure","formation_energy_per_atom"]
            )
    if len(docs) == 0:
          return None
    elif len(docs) == 1:
          pmg_structure = docs[0].structure
          ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure, msonable=False)
          return ase_atoms
    else:
          pmg_structure = min(docs,key=lambda x: x.formation_energy_per_atom).structure
          ase_atoms = AseAtomsAdaptor.get_atoms(pmg_structure, msonable=False)
          return ase_atoms
          #Should it be ase_atoms are some asemoolece

def random_3d_point_within_cell(v1, v2, v3):
        # Compute the volume of the parallelepiped
        volume = np.abs(np.dot(v1, np.cross(v2, v3)))
        
        # Generate random numbers
        r1, r2 = np.random.rand(2)
        r3 = np.random.rand() * (1 - r1 - r2)
        
        # Calculate the random point
        random_point = r1 * v1 + r2 * v2 + r3 * v3
        
        return random_point

def apply_modification(input_molecule,operation,mod_params):

    chemical_symbols = input_molecule.get_chemical_symbols()
    cell = input_molecule.get_cell()
    positions = input_molecule.get_positions()
    if operation == "substitute":
        atom1, atom2 = mod_params
        new_symbols = [atom2 if x == atom1 else x for x in chemical_symbols]
        new_positions = positions
    elif operation == "exchange":
        atom1, atom2 = mod_params
        new_symbols = [atom2 if x == atom1 else atom1 if x == atom2 else x for x in chemical_symbols]
        new_positions = positions
    elif operation == "add":
        atom = mod_params[0]
        new_symbols = chemical_symbols + [atom]
        new_positions = np.vstack(
                (
                    positions,
                    random_3d_point_within_cell(cell[0], cell[1], cell[2])
                )
            )
    elif operation == "remove":
        atom = mod_params[0]
        new_symbols = [x for x in chemical_symbols if x != atom]
        new_positions = positions[[x != atom for x in chemical_symbols]]

    else:
        ValueError(f"Invalid modification type: {operation}")

    new_molecule = Atoms(
                symbols=new_symbols,
                positions=new_positions,
                cell=cell,
                pbc=(True, True, True)
            )
    return new_molecule

    
def apply_random_mod(input_formula):
    modifications = ['substitute','exchange','add','remove']
    input_mol = convert_to_atom(input_formula)

    mod = random.choice(modifications)
    present_symbols = list(set(input_mol.get_chemical_symbols()))
    if mod == "add": 
        mod_param = [random.choice(chemical_symbols)]
    elif mod == "remove":
        mod_param = [random.choice(present_symbols)]
    elif mod == "exchange":
        mod_param = random.sample(present_symbols,2)
    else:
        mod_param =[random.choice(present_symbols),random.choice(chemical_symbols)]

    new_mol = apply_modification(input_mol,mod,mod_param)

    output_formula = new_mol.get_chemical_formula('metal')
    return (input_formula,output_formula,[mod] + mod_param)

if __name__ == "__main__":
     out = apply_random_mod('SrTiO3')
     print(out)      