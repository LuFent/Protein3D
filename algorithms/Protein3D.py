
from .Arcitecture import CIHBS
from .CustomPDBParser import PDBParser
from Bio import PDB as pdb
from Bio.PDB import Selection
import json
from itertools import chain, count
import numpy as np


class StructureVisualisation:
    cihbsObj = CIHBS.CIHBS()
    PDBParser = PDBParser(QUIET=True)

    def __init__(self, structure_id, structure_file):

        structure = self.PDBParser.get_structure(structure_id, structure_file)
        self.structure = structure

        self.atoms_ids = [atom.get_serial_number() for atom in
                          Selection.unfold_entities(self.structure, 'A')]

        self.cihbsObj.setParameters(structure=structure)

        #for atom in structure.get_atoms():
        #    print(f"{atom.get_serial_number()} {atom.get_name()}")

        self.extra_bonds = []
        self.hide_atoms = []
        self.show_atoms = []

    def flush_mask(self):
        self.hide_atoms = []
        self.extra_bonds = []
        self.show_atoms = []

    def get_mask(self):
        return {
            "hideAtoms": self.hide_atoms,
            "showAtoms": self.show_atoms,
            "extraBonds": self.extra_bonds,
        }

    def get_alpha_carbon_chain(self):
        self.flush_mask()

        show_ids = []
        extra_bonds = []

        previous_atom = None

        for model in self.structure:
            # Iterate through each chain in the model
            for chain in model:
                # Iterate through each residue in the chain
                for residue in chain:
                    # Iterate through each atom in the residue
                    for atom in residue:
                        if atom.get_name() == "CA":
                            show_ids.append(atom.get_serial_number())
                            if previous_atom is None:
                                previous_atom = atom.get_serial_number()
                            else:
                                extra_bonds.append((previous_atom, atom.get_serial_number()))
                                previous_atom = atom.get_serial_number()
                previous_atom = None

        self.show_atoms = show_ids
        self.extra_bonds = extra_bonds

    def getCIHBS(self):
        self.flush_mask()

        innerCIHBS = self.cihbsObj.checkInnerGroups(self.structure)
        physicalOperators = self.cihbsObj.connectPhysicalOperators()
        outerCIHBS = self.cihbsObj.connectResidueCIHBS()

        totalCIHBS = np.concatenate((physicalOperators, outerCIHBS), axis=0)
        print("totalCIHBS", totalCIHBS)

        CIHBS_ids = set(chain.from_iterable(innerCIHBS))
        self.show_atoms = list(CIHBS_ids)

    def test_alg(self):
        self.flush_mask()
        hide_atoms = []
        for model in self.structure:
            # Iterate through each chain in the model
            for chain in model:
                # Iterate through each residue in the chain
                for residue in chain:
                    if residue.get_resname() != "HOH":
                        continue
                    for atom in residue:
                        hide_atoms.append(atom.get_serial_number())

        self.hide_atoms = hide_atoms


    def execute_algorithm(self, alg):
        if alg == "alpha_carbon_skeleton":
            self.get_alpha_carbon_chain()
            return self.get_mask()

        elif alg == "CIHBS":
            self.getCIHBS()
            return self.get_mask()

        elif alg == "test_alg":
            self.test_alg()
            return self.get_mask()

        else:
            return None

